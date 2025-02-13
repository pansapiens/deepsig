#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import argparse
import json
from pathlib import Path
import multiprocessing as mp

DESC = "DeepSig: Predictor of signal peptides in proteins"

import deepsig
import deepsig.deepsigconfig as cfg
from deepsig import workenv
from deepsig.helpers import setUpTFCPU, readdata, printDate, write_gff_output, get_json_output, detectsp, predictsp, write_processed_sequences, write_noss_sequences

pclasses = {2: 'SignalPeptide', 1: 'Transmembrane', 0: 'Other'}

if('DEEPSIG_ROOT' in os.environ):
  try:
      deepsig_root = os.environ['DEEPSIG_ROOT']
      deepsig_root_path = Path(deepsig_root).resolve()
      if(not deepsig_root_path.is_dir()):
          raise IOError()
      elif(not deepsig_root_path.resolve(cfg.DNN_MODEL_DIR).is_dir()):
          raise IOError()
      elif(not deepsig_root_path.resolve(cfg.CRF_MODEL_DIR).is_dir()):
          raise IOError()
      else:
        sys.path.append(deepsig_root)
  except:
      sys.exit(f'ERROR: wrong DeepSig root path! DEEPSIG_ROOT={deepsig_root}')
else:
  sys.exit("ERROR: required environment variable 'DEEPSIG_ROOT' is not set")

def main():
  parser = argparse.ArgumentParser(description=DESC)
  parser.add_argument("-f", "--fasta",
                      help = "The input multi-FASTA file name",
                      dest = "fasta", required = True)
  parser.add_argument("-o", "--outf",
                      help = "The output file",
                      dest = "outf", required = True)
  parser.add_argument("-k", "--organism",
                      help = "The organism the sequences belongs to",
                      choices=['euk', 'gramp', 'gramn'],
                      dest = "organism", required = True)
  parser.add_argument("-m", "--outfmt",
                      help = "The output format: json or gff3 (default)",
                      choices=['json', 'gff3'], required = False, default = "gff3")

  parser.add_argument("--output-processed",
                      help = "Output FASTA file containing mature sequences with signal peptides removed",
                      dest = "processed_fasta",
                      required = False)

  parser.add_argument("--output-noss",
                      help = "Output FASTA file containing sequences without predicted signal peptides",
                      dest = "noss_fasta",
                      required = False)

  parser.add_argument('-t', '--threads', action='store', type=int, default=mp.cpu_count(), help='Number of threads to use (default = number of available CPUs)')

  parser.add_argument('--version', action='version', version=deepsig.__version__)

  ns = parser.parse_args()

  setUpTFCPU(ns.threads)

  protein_jsons = []
  try:
    we = workenv.TemporaryEnv()
    printDate("Reading input data")
    X, records = readdata(ns.fasta, cfg.NTERM)
    accs = [record.id for record in records]
    seqs = [str(record.seq) for record in records]
    printDate("Read %d protein sequences" % len(accs))
    printDate("Detecting signal peptides")
    Y, Ytm, Ync, cls, Ytm_norm, Ync_norm = detectsp(X, ns.organism)
    n_signal = cls.count(2)
    printDate("Detected %d signal peptides" % n_signal)
    printDate("Predicting cleavage sites")
    cleavage = predictsp(X, cls, ns.organism, we, cpu=ns.threads)
    printDate("Writing results to output file")
    ofs = open(ns.outf, 'w')
    for i in range(len(accs)):
      if cls[i] == 2:
        reliability = Y[i]
      elif cls[i] == 1:
        reliability = Ytm_norm[i]
      else:
        reliability = Ync_norm[i]
      if ns.outfmt == "gff3":
        write_gff_output(accs[i], seqs[i], ofs, pclasses[cls[i]], reliability, cleavage[i])
      else:
        acc_json = get_json_output(accs[i], seqs[i], pclasses[cls[i]], reliability, cleavage[i])
        protein_jsons.append(acc_json)
    ofs.close()
    if ns.outfmt == "json":
      ofs = open(ns.outf, 'w')
      json.dump(protein_jsons, ofs, indent=5)
      ofs.close()

    if ns.processed_fasta:
      printDate("Writing processed sequences to %s" % ns.processed_fasta)
      write_processed_sequences(records, cls, cleavage, ns.processed_fasta)

    if ns.noss_fasta:
      printDate("Writing %d sequences without signal peptides to %s" % (len(accs) - n_signal, ns.noss_fasta))
      write_noss_sequences(records, cls, ns.noss_fasta)
  except:
    printDate("Errors occured during execution")
    printDate("Leaving outdir unchanged")
    raise
  else:
    we.destroy()
    pass
  sys.exit(0)

if __name__ == "__main__":
  main()
