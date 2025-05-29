

import math
import json
from pathlib import Path
import logging
import re
logging.basicConfig(level=logging.DEBUG)
#logging.basicConfig(level=logging.INFO)


#!/usr/bin/env python
# vim: syntax=python tabstop=4 expandtab


import os
import sys
import re
import glob
import shutil
import gzip
import subprocess
from Bio import SwissProt
from Bio import SeqIO

class GTF:
    g_dic = dict()
    t_dic = dict()
    p_dic = dict()
    g2t_dic = dict()
    t2g_dic = dict()
    is_addVersion_dic = dict()

    def add_gene_biotype(self, gene_id, text):
        self.g_dic[gene_id]['biotype'] = text
    def add_gene_name(self, gene_id, text):
        self.g_dic[gene_id]['symbol'] = text
    def add_gene_version(self, gene_id, text):
        self.g_dic[gene_id]['version'] = text

    def add_tran_biotype(self, tran_id, text):
        self.t_dic[tran_id]['biotype'] = text
    def add_tran_name(self, tran_id, text):
        self.t_dic[tran_id]['symbol'] = text
    def add_tran_version(self, tran_id, text):
        self.t_dic[tran_id]['version'] = text

    def add_prot_biotype(self, prot_id, text):
        self.p_dic[prot_id]['biotype'] = text
    def add_prot_name(self, prot_id, text):
        self.p_dic[prot_id]['symbol'] = text
    def add_prot_version(self, prot_id, text):
        self.p_dic[prot_id]['version'] = text

    def add_addVersionId(self, _id, _version):
        self.is_addVersion_dic.setdefault('{0}.{1}'.format(_id,_version), _id)

    def init_dics_from_gtf(self, gtf_fn):
        for line in open(gtf_fn):
            if line.startswith('#'):
                continue
            items = line.rstrip('\n').split('\t')
            attr_dic = dict()
            attr_s = [x.strip() for x in items[8].split(';') if x]
            for attr in attr_s:
                key = attr.split()[0]
                value_s = re.findall(r'\s\"(.*?)\"', attr)
                if value_s:
                    value = value_s[0]
                else:
                    continue
                attr_dic.setdefault(key, value)

            if 'gene_id' in attr_dic:
                gene_id = attr_dic['gene_id']
                self.g_dic.setdefault(gene_id, {}).setdefault('biotype', 'NA')
                self.g_dic.setdefault(gene_id, {}).setdefault('symbol', 'NA')
                self.g_dic.setdefault(gene_id, {}).setdefault('description', 'NA')
                self.g_dic.setdefault(gene_id, {}).setdefault('version', 'NA')
                if 'gene_biotype' in attr_dic:
                    self.add_gene_biotype(gene_id, attr_dic['gene_biotype'])
                if 'gene_name' in attr_dic:
                    self.add_gene_name(gene_id, attr_dic['gene_name'])
                if 'gene_version' in attr_dic:
                    self.add_gene_version(gene_id, attr_dic['gene_version'])
                    self.add_addVersionId(gene_id, attr_dic['gene_version'])

                self.g2t_dic.setdefault(gene_id, [])

            if 'transcript_id' in attr_dic:
                tran_id = attr_dic['transcript_id']
                self.t_dic.setdefault(tran_id, {}).setdefault('biotype', 'NA')
                self.t_dic.setdefault(tran_id, {}).setdefault('symbol', 'NA')
                self.t_dic.setdefault(tran_id, {}).setdefault('description', 'NA')
                self.t_dic.setdefault(tran_id, {}).setdefault('version', 'NA')
                if 'transcript_biotype' in attr_dic:
                    self.add_tran_biotype(tran_id, attr_dic['transcript_biotype'])
                if 'transcript_name' in attr_dic:
                    self.add_tran_name(tran_id, attr_dic['transcript_name'])
                if 'transcript_version' in attr_dic:
                    self.add_tran_version(tran_id, attr_dic['transcript_version'])
                    self.add_addVersionId(tran_id, attr_dic['transcript_version'])

                if not tran_id in self.g2t_dic[gene_id]:
                    self.g2t_dic[gene_id].append(tran_id)

                self.t2g_dic.setdefault(tran_id, gene_id)

            if 'protein_id' in attr_dic:
                prot_id = attr_dic['protein_id']
                self.p_dic.setdefault(prot_id, {}).setdefault('biotype', 'NA')
                self.p_dic.setdefault(prot_id, {}).setdefault('symbol', 'NA')
                self.p_dic.setdefault(prot_id, {}).setdefault('description', 'NA')
                self.p_dic.setdefault(prot_id, {}).setdefault('version', 'NA')
                if 'protein_biotype' in attr_dic:
                    self.add_prot_biotype(prot_id, attr_dic['protein_biotype'])
                if 'protein_name' in attr_dic:
                    self.add_prot_name(prot_id, attr_dic['protein_name'])
                if 'protein_version' in attr_dic:
                    self.add_prot_version(prot_id, attr_dic['protein_version'])
                    self.add_addVersionId(prot_id, attr_dic['protein_version'])


class GetAnnoUniprot(GTF):
    def __init__(self, wkdir, gtf, tsv, dat):
        self.init_dics_from_gtf(gtf)

        self.wkdir = wkdir
        self.tsv = tsv
        self.dat = dat

        self.xref_dic = dict()
        self.id_dic = dict()
        self.parse_tsv()
        self.parse_dat()

        self.write_anno('gene', 'go')
        self.write_anno('transcript', 'go')
        self.write_anno('gene', 'ipr')
        self.write_anno('transcript', 'ipr')
        self.write_anno('gene', 'pfam')
        self.write_anno('transcript', 'pfam')

    def parse_tsv(self):
        for line in open(self.tsv):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['gene_stable_id']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            gene_id = items[idx_dic['gene_stable_id']]
            if not gene_id in self.g_dic:
                continue
            tran_id = items[idx_dic['transcript_stable_id']]
            if not tran_id in self.t_dic:
                continue

            xref = items[idx_dic['xref']]

            db_name = items[idx_dic['db_name']]
            if db_name in ['Uniprot/SWISSPROT', 'Uniprot/SPTREMBL']:
                self.xref_dic.setdefault(xref, {}).setdefault('go',{})
                self.xref_dic.setdefault(xref, {}).setdefault('ipr',{})
                self.xref_dic.setdefault(xref, {}).setdefault('pfam',{})

                self.id_dic.setdefault('gene',{}).setdefault(gene_id, {}).setdefault(xref, {})
                self.id_dic.setdefault('transcript',{}).setdefault(tran_id, {}).setdefault(xref, {})
            #

    def parse_dat(self):
        for record in SwissProt.parse(open(self.dat)):
            xref = record.accessions[0]
            if not xref in self.xref_dic:
                continue
            for units in record.cross_references:
                if units[0] in ['GO']:
                    go_id = units[1]
                    go_desc = units[2]
                    go_evidence = units[3]
                    self.xref_dic[xref]['go'].setdefault(go_id, go_desc)
                elif units[0] in ['InterPro']:
                    ipr_id = units[1]
                    ipr_desc = units[2]
                    self.xref_dic[xref]['ipr'].setdefault(ipr_id, ipr_desc)
                elif units[0] in ['Pfam']:
                    pfam_id = units[1]
                    pfam_desc = units[2]
                    self.xref_dic[xref]['pfam'].setdefault(pfam_id, pfam_desc)
                else:
                    pass
            #

    def write_anno(self, level, db):
        outfn = os.path.join(self.wkdir, 'anno.{0}2{1}.tsv'.format(level, db))
        outfh = open(outfn, 'w')
        outfh.write('{0}\n'.format('\t'.join(['ID','ANNO_ID','ANNO_DESC'])))
        for level_id, xref_dic in sorted(self.id_dic[level].items()):
            for xref_id in xref_dic:
                anno_dic = self.xref_dic[xref_id][db]
                for anno_id, anno_desc in sorted(anno_dic.items()):
                    items = [level_id, anno_id, anno_desc]
                    outfh.write('{0}\n'.format('\t'.join(items)))
        outfh.close()


class GetDesc(GTF):
    def __init__(self, wkdir, gtf, mysqlgene, mysqltran):
        self.init_dics_from_gtf(gtf)
        self.add_desc(mysqlgene)
        #self.add_gene_desc(mysqlgene)
        #self.add_tran_desc(mysqltran)

        self.writer(self.g_dic, os.path.join(wkdir, 'anno.gene2desc.tsv'))
        self.writer(self.t_dic, os.path.join(wkdir, 'anno.transcript2desc.tsv'))

    def add_desc(self, fn):
        for line in open(fn):
            items = line.rstrip('\n').split('\t')
            desc = items[9]
            g_id = items[12]
            if not g_id in self.g_dic:
                print('{0} is not in gtf'.format(g_id))
                continue
            self.g_dic[g_id]['description'] = desc
            for t_id in self.g2t_dic[g_id]:
                self.t_dic[t_id]['description'] = desc

    def writer(self, _dic, outfn):
        outfh = open(outfn, 'w')
        outfh.write('{0}\n'.format('\t'.join(['ID','BIOTYPE','SYMBOL','DESCRIPTION'])))
        for _id, info_dic in sorted(_dic.items()):
            items = [_id]
            items.append(info_dic['biotype'])
            items.append(info_dic['symbol'])
            items.append(info_dic['description'])
            outfh.write('{0}\n'.format('\t'.join(items)))
        outfh.close()

    def _find_pattern(self, pattern, string):
        m = re.findall(pattern, string)
        if m:
            return m[0]
        else:
            return ''

    def _parse_fasta(self):
        for record in SeqIO.parse(open(self.fasta), 'fasta'):
            transcript = record.id
            _desc = record.description
            chromosome = self.find_pattern(r'\schromosome\:(.*?)\sgene\:', _desc)
            gene = self.find_pattern(r'\sgene\:(.*?)\sgene_biotype:', _desc)
            biotype = self.find_pattern(r'\sgene_biotype\:(.*?)\stranscript_biotype:', _desc)
            symbol = self.find_pattern(r'\sgene_symbol\:(.*?)\sdescription:', _desc)
            description = self.find_pattern(r'\sdescription\:(.*?)$', _desc)

            if gene in self.id_version_dic:
                gene = self.id_version_dic[gene]
            if transcript in self.id_version_dic:
                transcript = self.id_version_dic[transcript]

            if gene in self.g_dic:
                self.g_dic[gene]['biotype'] = biotype
                self.g_dic[gene]['symbol'] = symbol
                self.g_dic[gene]['description'] = description

            if transcript in self.t_dic:
                self.t_dic[transcript]['biotype'] = biotype
                self.t_dic[transcript]['symbol'] = symbol
                self.t_dic[transcript]['description'] = description

class GetMapG2TGTF(GTF):
    def __init__(self, wkdir, gtf):
        self.init_dics_from_gtf(gtf)
        self.writer(os.path.join(wkdir, 'map.gene2transcript.tsv'))

    def writer(self, outfn):
        outfh = open(outfn, 'w')
        for gene_id, tran_id_s in sorted(self.g2t_dic.items()):
            if len(tran_id_s) in [0]:
                print('{0} has not transcripts'.format(gene_id))
                sys.exit()
            for tran_id in sorted(tran_id_s):
                g_biotype = self.g_dic[gene_id]['biotype']
                t_biotype = self.t_dic[tran_id]['biotype']
                outfh.write('{0}\n'.format('\t'.join([gene_id, g_biotype, tran_id, t_biotype])))
        outfh.close()

class CleanupFasta(GTF):
    def __init__(self, wkdir, target, gtf, infasta, outfasta):
        self.init_dics_from_gtf(gtf)

        self.makefasta(target, infasta, outfasta)

    def makefasta(self, target, infasta, outfasta):
        outfh = open(outfasta, 'w')
        for record in SeqIO.parse(open(infasta), 'fasta'):
            if record.id in self.is_addVersion_dic:
                record.id = self.is_addVersion_dic[record.id]
            if target in ['cds', 'cdna']:
                if record.id in self.t_dic:
                    SeqIO.write(record, outfh, 'fasta')
            elif target in ['pep']:
                if record.id in self.p_dic:
                    SeqIO.write(record, outfh, 'fasta')
            else:
                pass
        outfh.close()


class Fastp:
    def __init__(self, wkdir, outprefix):
        self.wkdir = Path(wkdir)
        self.outprefix = outprefix
        self.json_path_dic = dict()

    def find_json_path(self, target_dir):
        for json_path in Path(target_dir).rglob("*.fastp.json"):
            sample_name = json_path.parts[-1].replace(".fastp.json", "")
            logging.info(f"Found {str(json_path)} of {sample_name}")
            self.json_path_dic.setdefault(sample_name, json.load(json_path.open("r")))
        return 1

    def write_summary_tsv(self, tag):
        outfp = self.wkdir / f"{self.outprefix}.xls"
        logging.info(f"writing fastp summary table... {str(outfp)}")
        outfh = outfp.open("w")
        headers = ["Samples"]
        headers.append("ReadsCount")
        headers.append("BasesCount")
        headers.append("BasesCountGb")
        headers.append("GCRate")
        headers.append("Q20BaseRate")
        headers.append("Q30BaseRate")
        headers.append("DupRate")
        headers.append("InsertSize")
        headers.append("R1MeanLen")
        headers.append("R2MeanLen")
        headers.append("GoodReadRate")
        outfh.write('{0}\n'.format("\t".join(headers)))
        for sample_name, info_dic in sorted(self.json_path_dic.items()):
            items = [sample_name]
            total_reads = info_dic["summary"][tag]["total_reads"]
            items.append(total_reads)
            items.append(info_dic["summary"][tag]["total_bases"])
            items.append(round(int(info_dic["summary"][tag]["total_bases"])*0.000000001, 3))
            items.append(round(float(info_dic["summary"][tag]["gc_content"])*100, 3))
            items.append(round(float(info_dic["summary"][tag]["q20_rate"])*100, 3))
            items.append(round(float(info_dic["summary"][tag]["q30_rate"])*100, 3))
            items.append(round(float(info_dic["duplication"]["rate"])*100, 3))
            if "insert_size" in info_dic:
                items.append(info_dic["insert_size"]["peak"])
            else:
                items.append("-")
            items.append(info_dic["summary"][tag]["read1_mean_length"])
            if "read2_mean_length" in info_dic["summary"][tag]:
                items.append(info_dic["summary"][tag]["read2_mean_length"])
            else:
                items.append("-")
            items.append(round(int(info_dic["filtering_result"]["passed_filter_reads"])/float(total_reads)*100, 3))
            outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
        outfh.close()
        return 1


class Tophat:
    def __init__(self, wkdir, outprefix, target_dir):
        self.wkdir = Path(wkdir)
        self.input_read_type = None

        # parse align_summary.txt files
        self.align_summary_path_dic = dict()
        self.align_summary_info_dic = dict()
        self.find_align_summary_path(target_dir)
        self.parse_align_summary_path()
        self.write_align_summary_tsv(Path(f"{outprefix}.xls"))

    def write_align_summary_tsv(self, outfp):
        logging.info(f"writing align summary table... {str(outfp)}")
        outfh = outfp.open("w")
        headers = ["SampleName"]
        headers.append("InputCount")
        headers.append("MappedCount")
        headers.append("MappedRate")
        headers.append("MultiMapCount")
        headers.append("MultiMapRate")
        headers.append("AlignedAsPairCount")
        headers.append("AlignedAsPairRate")
        headers.append("AlignedAsPair-MultiRate")
        headers.append("AlignedAsPair-DiscordantRate")
        headers.append("AlignedAsPair-ConcordantRate")
        outfh.write("{0}\n".format('\t'.join([str(x) for x in headers])))
        for sample_name, info_dic in self.align_summary_info_dic.items():
            items = [sample_name]
            items.append(info_dic["total"]["input"])
            items.append(info_dic["total"]["mapped"])
            items.append(round((info_dic["total"]["mapped"]/(info_dic["total"]["input"]*1.0))*100, 2))
            items.append(info_dic["total"]["multimap"])
            items.append(round((info_dic["total"]["multimap"]/(info_dic["total"]["mapped"]*1.0))*100, 2))
            items.append(info_dic["total"]["aligned-pair"])
            if self.input_read_type in ["paired-end"]:
                items.append(round((info_dic["total"]["aligned-pair"]/(info_dic["total"]["mapped"]*1.0))*100, 2))
            else:
                items.append(0)
            if self.input_read_type in ["paired-end"]:
                items.append(round((info_dic["total"]["aligned-pair-multi"]/(info_dic["total"]["aligned-pair"]*1.0))*100, 2))
            else:
                items.append(0)
            if self.input_read_type in ["paired-end"]:
                items.append(round((info_dic["total"]["aligned-pair-discordant"]/(info_dic["total"]["aligned-pair"]*1.0))*100, 2))
            else:
                items.append(0)
            if self.input_read_type in ["paired-end"]:
                items.append(round(info_dic["total"]["aligned-pair-concordant-rate"], 2))
            else:
                items.append(0)
            logging.info(items)
            outfh.write("{0}\n".format('\t'.join([str(x) for x in items])))
        outfh.close()
        return 1


    def parse_align_summary_path(self):
        for sample_name, align_summary_path in self.align_summary_path_dic.items():
            logging.debug(sample_name)
            metrix = list()
            for line in align_summary_path.open():
                items = line.strip().split()
                if not items:
                    continue
                metrix.append(items)
            logging.debug(len(metrix))
            for items in metrix:
                logging.debug(items)

            _dic = dict()
            _dic.setdefault("r1", {}).setdefault("input", 0)
            _dic.setdefault("r1", {}).setdefault("mapped", 0)
            _dic.setdefault("r1", {}).setdefault("multimap", 0)
            _dic.setdefault("r2", {}).setdefault("input", 0)
            _dic.setdefault("r2", {}).setdefault("mapped", 0)
            _dic.setdefault("r2", {}).setdefault("multimap", 0)
            _dic.setdefault("total", {}).setdefault("input", 0)
            _dic.setdefault("total", {}).setdefault("mapped", 0)
            _dic.setdefault("total", {}).setdefault("multimap", 0)
            _dic.setdefault("total", {}).setdefault("aligned-pair", 0)
            _dic.setdefault("total", {}).setdefault("aligned-pair-multi", 0)
            _dic.setdefault("total", {}).setdefault("aligned-pair-discordant", 0)
            _dic.setdefault("total", {}).setdefault("aligned-pair-concordant-rate", 0.0)
            if len(metrix) in [13]: # align_summary for paired-end reads
                self.input_read_type = "paired-end"

                _dic["r1"]["input"] = int(metrix[1][2])
                _dic["r1"]["mapped"] = int(metrix[2][2])
                _dic["r1"]["multimap"] = int(metrix[3][2])
                _dic["r2"]["input"] = int(metrix[5][2])
                _dic["r2"]["mapped"] = int(metrix[6][2])
                _dic["r2"]["multimap"] = int(metrix[7][2])
                _dic["total"]["input"] = _dic["r1"]["input"] + _dic["r2"]["input"]
                _dic["total"]["mapped"] = _dic["r1"]["mapped"] + _dic["r2"]["mapped"]
                _dic["total"]["multimap"] = _dic["r1"]["multimap"] + _dic["r2"]["multimap"]
                _dic["total"]["aligned-pair"] = int(metrix[9][2])*2
                _dic["total"]["aligned-pair-multi"] = int(metrix[10][2])*2
                _dic["total"]["aligned-pair-discordant"] = int(metrix[11][0])*2
                _dic["total"]["aligned-pair-concordant-rate"] = float(metrix[12][0].replace("%", ""))
            elif len(metrix) in [5]: # align_summary for single-end reads
                self.input_read_type = "single-end"

                _dic["r1"]["input"] = int(metrix[1][2])
                _dic["r1"]["mapped"] = int(metrix[2][2])
                _dic["r1"]["multimap"] = int(metrix[3][2])
                _dic["total"]["input"] = _dic["r1"]["input"]
                _dic["total"]["mapped"] = _dic["r1"]["mapped"]
                _dic["total"]["multimap"] = _dic["r1"]["multimap"]
            else:
                logging.critical(f"len(metrix) is {len(metrix)}")
                logging.critical("It is need to dev for this.")

            self.align_summary_info_dic.setdefault(sample_name, _dic)
        logging.debug(self.align_summary_info_dic)
        return 1

    def find_align_summary_path(self, tophat_dir):
        for align_summary_path in Path(tophat_dir).rglob("*/align_summary.txt"):
            sample_name = align_summary_path.parts[-2]
            logging.info(f"Found {str(align_summary_path)} of {sample_name}")
            self.align_summary_path_dic.setdefault(sample_name, align_summary_path)
        return 1


class Cufflinks:
    def __init__(self, wkdir, outprefix, gene2desc, tran2desc, gene2go, tran2go):
        self.wkdir = Path(wkdir)
        self.sample_names = list()
        self.gene2desc_fp = Path(gene2desc)
        self.tran2desc_fp = Path(tran2desc)
        self.gene2go_fp = Path(gene2go)
        self.tran2go_fp = Path(tran2go)

        # setting vars for annotation
        self.header_names = list()
        self.anno_names = list()
        _names = [("locus", "Locus"),
                  ("name", "GeneName"),
                  ("symbol", "GeneSymbol"),
                  ("description", "Description"),
                  ("biotype", "Biotype")]
        for anno_name, header_name in _names:
            self.anno_names.append(anno_name)
            self.header_names.append(header_name)

        # parse gene.fpkm_tracking files | cufflinks with -G option
        self.gene_path_dic = dict()
        self.gene_info_dic = dict()
        self.gene_fpkm_dic = dict()
        self.gene_stat_dic = dict()


        # parse isoforms.fpkm_tracking files | cufflinks with -G option
        self.tran_path_dic = dict()
        self.tran_info_dic = dict()
        self.tran_fpkm_dic = dict()
        self.tran_stat_dic = dict()

    def write_stat_coexp(self, molecule_type, outfp):
        if molecule_type in ["genes"]:
            _stat_dic = self.gene_stat_dic
            _fpkm_dic = self.gene_fpkm_dic
            _info_dic = self.gene_info_dic
        elif molecule_type in ["isoforms"]:
            _stat_dic = self.tran_stat_dic
            _fpkm_dic = self.tran_fpkm_dic
            _info_dic = self.tran_info_dic

        logging.info(f"writing STAT of Co-Expressed gene table... {str(outfp)}")

        _stat_col_s = ["fpkm>=0.3", "fpkm>=1", "fpkm>=10", "fpkm>=100", "fpkm>=1000"]

        _coexp_dic = dict()
        for i in range(len(self.sample_names)+1):
            for _stat_col in _stat_col_s:
                _coexp_dic.setdefault(i, {}).setdefault(_stat_col, 0)

        for _id, info_dic in _fpkm_dic.items():

            if _info_dic[_id]['biotype'] not in ['protein_coding']:
                continue # counting only protein_coding

            _fpkm_s = list()
            for sample_name in self.sample_names:
                _fpkm_s.append(float(info_dic[sample_name]))
            for _stat_col in _stat_col_s:
                _fpkm_cutoff = 0.0
                if _stat_col in ["fpkm>=0.3"]:
                    _fpkm_cutoff = 0.3
                elif _stat_col in ["fpkm>=1"]:
                    _fpkm_cutoff = 1
                elif _stat_col in ["fpkm>=10"]:
                    _fpkm_cutoff = 10
                elif _stat_col in ["fpkm>=100"]:
                    _fpkm_cutoff = 100
                elif _stat_col in ["fpkm>=1000"]:
                    _fpkm_cutoff = 1000
                _isexp_s = ["O" if x >= _fpkm_cutoff else "X" for x in _fpkm_s]
                _coexp_dic[_isexp_s.count("O")][_stat_col] += 1

        headers = ["CoExpSampleCount"]
        headers.extend(_stat_col_s)
        outfh = outfp.open("w")
        outfh.write("{0}\n".format("\t".join([str(x) for x in headers])))
        for coexp_sample_count, info_dic in sorted(_coexp_dic.items(), reverse=True):
            items = [coexp_sample_count]
            for _stat_col in _stat_col_s:
                items.append(info_dic[_stat_col])
            outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
        outfh.close()
        return 1

    def write_stat_exp(self, molecule_type, outfp):
        if molecule_type in ["genes"]:
            _stat_dic = self.gene_stat_dic
        elif molecule_type in ["isoforms"]:
            _stat_dic = self.tran_stat_dic

        logging.info(f"writing STAT of Expressed gene table... {str(outfp)}")

        _stat_col_s = ["fpkm>=0.3", "fpkm>=1", "fpkm>=10",
                       "fpkm>=100", "fpkm>=1000"]
        headers = ["SampleNames"]
        headers.extend(_stat_col_s)
        outfh = outfp.open("w")
        outfh.write("{0}\n".format("\t".join([str(x) for x in headers])))
        for sample_name in self.sample_names:
            info_dic = _stat_dic[sample_name]
            items = [sample_name]
            for _stat_col in _stat_col_s:
                items.append(len(info_dic[_stat_col]))
            outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
        outfh.close()
        return 1


    def write_fpkm_tsv(self, molecule_type, fpkm_criteria, outfp):
        if molecule_type in ["genes"]:
            _info_dic = self.gene_info_dic
            _fpkm_dic = self.gene_fpkm_dic
        elif molecule_type in ["isoforms"]:
            _info_dic = self.tran_info_dic
            _fpkm_dic = self.tran_fpkm_dic

        logging.info(f"writing FPKM table... {str(outfp)}")

        headers = ["GeneID"]
        headers.extend(self.sample_names)
        headers.extend(self.header_names)

        outfh = outfp.open("w")
        outfh.write("{0}\n".format("\t".join([str(x) for x in headers])))
        for _id, _sub_dic in _fpkm_dic.items():
            items = [_id]
            fpkm_s = list()
            for sample_name in self.sample_names:
                _fpkm = _sub_dic[sample_name]
                fpkm_s.append(_fpkm)
                items.append(round(_fpkm, 4))
                #items.append(int(_fpkm))

            for anno_name in self.anno_names:
                if anno_name in _info_dic[_id]:
                    items.append(_info_dic[_id][anno_name])
                else:
                    items.append("None")

            if max(fpkm_s) >= fpkm_criteria:
                outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
        outfh.close()
        return 1

    def anno_gene2desc(self, colname, key):
        for line in self.gene2desc_fp.open():
            items = line.rstrip('\n').split('\t')
            if items[0] in ["ID"]:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            _id = items[idx_dic["ID"]]
            if _id in self.gene_info_dic:
                self.gene_info_dic[_id].setdefault(key, items[idx_dic[colname]])
        return 1

    def anno_xxxx2desc(self, molecule_type, xxxx2desc_fp, colname, key):
        if molecule_type in ["genes"]:
            _info_dic = self.gene_info_dic
        elif molecule_type in ["isoforms"]:
            _info_dic = self.tran_info_dic

        for line in xxxx2desc_fp.open():
            items = line.rstrip('\n').split('\t')
            if items[0] in ["ID"]:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            _id = items[idx_dic["ID"]]
            _info_dic.setdefault(_id, {}).setdefault(key, items[idx_dic[colname]])

        if molecule_type in ["genes"]:
            self.gene_info_dic = _info_dic
        elif molecule_type in ["isoforms"]:
            self.tran_info_dic = _info_dic

        return 1

    def anno_xxxx2go(self, molecule_type, xxxx2go_fp):
        if molecule_type in ["genes"]:
            _info_dic = self.gene_info_dic
        elif molecule_type in ["isoforms"]:
            _info_dic = self.tran_info_dic

        # add GO to headers
        for go_category in ["CC", "MF", "BP"]:
            if f"GO:{go_category}" not in self.header_names:
                self.header_names.append(f"GO:{go_category}")
            if f"GO:{go_category}" not in self.anno_names:
                self.anno_names.append(f"GO:{go_category}")

        for line in xxxx2go_fp.open():
            items = line.rstrip('\n').split('\t')
            if items[0] in ["ID"]:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            _id = items[idx_dic["ID"]]
            go_id = items[idx_dic["ANNO_ID"]]
            _go_desc = items[idx_dic["ANNO_DESC"]]
            _go_category = _go_desc.split(':')[0]
            go_category = None
            if _go_category in ["C"]:
                go_category = "CC"
            elif _go_category in ["F"]:
                go_category = "MF"
            elif _go_category in ["P"]:
                go_category = "BP"
            go_term = _go_desc.split(':')[1]

            # init GO dic
            if _id in _info_dic:
                _info_dic[_id].setdefault("GO:CC", {})
                _info_dic[_id].setdefault("GO:MF", {})
                _info_dic[_id].setdefault("GO:BP", {})
                _info_dic[_id][f"GO:{go_category}"].setdefault(go_id, go_term)
            else:
                logging.warning(f"WARNING: {_id} is not in {molecule_type} info dic")


        if molecule_type in ["genes"]:
            self.gene_info_dic = _info_dic
        elif molecule_type in ["isoforms"]:
            self.tran_info_dic = _info_dic

        return 1

    def parse_fpkm_tracking_path(self, molecule_type):
        if molecule_type in ["genes"]:
            _path_dic = self.gene_path_dic
            _info_dic = self.gene_info_dic
        elif molecule_type in ["isoforms"]:
            _path_dic = self.tran_path_dic
            _info_dic = self.tran_info_dic

        _fpkm_dic = dict()
        _stat_dic = dict()
        for sample_name, fpkm_tracking_path in _path_dic.items():
            for line in fpkm_tracking_path.open():
                items = line.rstrip('\n').split('\t')
                if items[0] in ["tracking_id"]:
                    idx_dic = dict()
                    for idx, item in enumerate(items):
                        idx_dic.setdefault(item, idx)
                    continue
                _id = items[idx_dic["tracking_id"]]
                _name = items[idx_dic["gene_short_name"]]
                _locus = items[idx_dic["locus"]]
                _fpkm = float(items[idx_dic["FPKM"]])

                _info_dic.setdefault(_id, {}).setdefault("name", _name)
                _info_dic.setdefault(_id, {}).setdefault("locus", _locus)
                _fpkm_dic.setdefault(_id, {}).setdefault(sample_name, _fpkm)
                _stat_dic.setdefault(sample_name, {}).setdefault("fpkm>=0.3", [])
                _stat_dic.setdefault(sample_name, {}).setdefault("fpkm>=1", [])
                _stat_dic.setdefault(sample_name, {}).setdefault("fpkm>=10", [])
                _stat_dic.setdefault(sample_name, {}).setdefault("fpkm>=100", [])
                _stat_dic.setdefault(sample_name, {}).setdefault("fpkm>=1000", [])

                if float(_fpkm) >= 1000 and _info_dic[_id]['biotype'] in ['protein_coding']: # counting only protein_coding
                    _stat_dic[sample_name]["fpkm>=1000"].append(_id)
                if float(_fpkm) >= 100 and _info_dic[_id]['biotype'] in ['protein_coding']: # counting only protein_coding
                    _stat_dic[sample_name]["fpkm>=100"].append(_id)
                if float(_fpkm) >= 10 and _info_dic[_id]['biotype'] in ['protein_coding']: # counting only protein_coding
                    _stat_dic[sample_name]["fpkm>=10"].append(_id)
                if float(_fpkm) >= 1.0 and _info_dic[_id]['biotype'] in ['protein_coding']: # counting only protein_coding
                    _stat_dic[sample_name]["fpkm>=1"].append(_id)
                if float(_fpkm) >= 0.3 and _info_dic[_id]['biotype'] in ['protein_coding']: # counting only protein_coding
                    _stat_dic[sample_name]["fpkm>=0.3"].append(_id)


        if molecule_type in ["genes"]:
            self.gene_info_dic = _info_dic
            self.gene_fpkm_dic = _fpkm_dic
            self.gene_stat_dic = _stat_dic
        elif molecule_type in ["isoforms"]:
            self.tran_info_dic = _info_dic
            self.tran_fpkm_dic = _fpkm_dic
            self.tran_stat_dic = _stat_dic

        return 1

    def find_fpkm_tracking_path(self, molecule_type, cufflinks_dir):
        for fpkm_tracking_path in Path(cufflinks_dir).glob(f"*/{molecule_type}.fpkm_tracking"):
            sample_name = fpkm_tracking_path.parts[-2]
            if sample_name not in self.sample_names:
                self.sample_names.append(sample_name)
            logging.info(f"Found {str(fpkm_tracking_path)} of {sample_name}")

            if molecule_type in ["genes"]:
                self.gene_path_dic.setdefault(sample_name, fpkm_tracking_path)
            elif molecule_type in ["isoforms"]:
                self.tran_path_dic.setdefault(sample_name, fpkm_tracking_path)
        self.sample_names.sort()

        return 1


class Cuffnorm(Cufflinks):

    def get_sample_names_from_fpkm_table(self, fp):

        for line in fp.open():
            items = line.rstrip('\n').split('\t')
            if items[0] in ["tracking_id"]:
                _sample_names = items[1:]
            else:
                break

        for _sample_name in _sample_names:
            sample_name = re.sub("_0$", "", _sample_name)
            if sample_name not in self.sample_names:
                self.sample_names.append(sample_name)

        self.sample_names.sort()

        return 1

    def find_fpkm_table_path(self, molecule_type, cuffnorm_dir):

        fpkm_table_path = Path(cuffnorm_dir) / f"{molecule_type}.fpkm_table"
        self.get_sample_names_from_fpkm_table(fpkm_table_path)
        logging.info(f"Found samples {self.sample_names} in {str(fpkm_table_path)}")

        if molecule_type in ["genes"]:
            for sample_name in self.sample_names:
                self.gene_path_dic.setdefault(sample_name, fpkm_table_path)
        elif molecule_type in ["isoforms"]:
            for sample_name in self.sample_names:
                self.tran_path_dic.setdefault(sample_name, fpkm_table_path)

        return 1

    def parse_fpkm_table_path(self, molecule_type):
        if molecule_type in ["genes"]:
            _path_dic = self.gene_path_dic
            _info_dic = self.gene_info_dic
        elif molecule_type in ["isoforms"]:
            _path_dic = self.tran_path_dic
            _info_dic = self.tran_info_dic

        _fpkm_dic = dict()
        _stat_dic = dict()
        for sample_name, fpkm_table_path in _path_dic.items():
            for line in fpkm_table_path.open():
                items = line.rstrip('\n').split('\t')
                if items[0] in ["tracking_id"]:
                    idx_dic = dict()
                    for idx, item in enumerate(items):
                        idx_dic.setdefault(item, idx)
                    continue
                _id = items[idx_dic["tracking_id"]]
                _name = "-"
                _locus = "-"
                _fpkm = float(items[idx_dic[f"{sample_name}_0"]])

                _info_dic.setdefault(_id, {}).setdefault("name", _name)
                _info_dic.setdefault(_id, {}).setdefault("locus", _locus)
                _fpkm_dic.setdefault(_id, {}).setdefault(sample_name, _fpkm)
                _stat_dic.setdefault(sample_name, {}).setdefault("fpkm>=0.3", [])
                _stat_dic.setdefault(sample_name, {}).setdefault("fpkm>=1", [])
                _stat_dic.setdefault(sample_name, {}).setdefault("fpkm>=10", [])
                _stat_dic.setdefault(sample_name, {}).setdefault("fpkm>=100", [])
                _stat_dic.setdefault(sample_name, {}).setdefault("fpkm>=1000", [])

                if float(_fpkm) >= 1000 and _info_dic[_id]['biotype'] in ['protein_coding']: # counting only protein_coding
                    _stat_dic[sample_name]["fpkm>=1000"].append(_id)
                if float(_fpkm) >= 100 and _info_dic[_id]['biotype'] in ['protein_coding']: # counting only protein_coding
                    _stat_dic[sample_name]["fpkm>=100"].append(_id)
                if float(_fpkm) >= 10 and _info_dic[_id]['biotype'] in ['protein_coding']: # counting only protein_coding
                    _stat_dic[sample_name]["fpkm>=10"].append(_id)
                if float(_fpkm) >= 1.0 and _info_dic[_id]['biotype'] in ['protein_coding']: # counting only protein_coding
                    _stat_dic[sample_name]["fpkm>=1"].append(_id)
                if float(_fpkm) >= 0.3 and _info_dic[_id]['biotype'] in ['protein_coding']: # counting only protein_coding
                    _stat_dic[sample_name]["fpkm>=0.3"].append(_id)


        if molecule_type in ["genes"]:
            self.gene_info_dic = _info_dic
            self.gene_fpkm_dic = _fpkm_dic
            self.gene_stat_dic = _stat_dic
        elif molecule_type in ["isoforms"]:
            self.tran_info_dic = _info_dic
            self.tran_fpkm_dic = _fpkm_dic
            self.tran_stat_dic = _stat_dic

        return 1


class Cuffdiff:
    def __init__(self, wkdir, outprefix, log2fc_criteria, p_criteria, q_criteria):
        self.wkdir = Path(wkdir)
        self.outprefix = outprefix
        self.comp_names = list()
        self._diff_col_s = ['Status', 'CtrlName','CaseName',
                            'CtrlFPKM','CaseFPKM',
                            'Log2FC','p','q', 'UpDown','IsDEG']
        self.gene_info_col_s = ['GeneName', 'Locus']
        self.tran_info_col_s = ['GeneName', 'Locus']
        self.log2fc_criteria = log2fc_criteria
        self.p_criteria = p_criteria
        self.q_criteria = q_criteria

        # parse gene_exp.diff files
        self.gene_path_dic = dict()
        self.gene_info_dic = dict()
        self.gene_diff_dic = dict()
        self.gene_stat_dic = dict()
        self.find_exp_diff_path("gene", "anl/cuffdiff")
        self.parse_exp_diff_path("gene")
        self.stat_exp_diff("gene")

        # parse isoform_exp_diff files
        self.tran_path_dic = dict()
        self.tran_info_dic = dict()
        self.tran_diff_dic = dict()
        self.find_exp_diff_path("isoform", "anl/cuffdiff")
        self.parse_exp_diff_path("isoform")
        self.stat_exp_diff("isoform")

    def find_exp_diff_path(self, molecule_type, cuffdiff_dir):
        for exp_diff_path in Path(cuffdiff_dir).glob(f"*/{molecule_type}_exp.diff"):
            comp_name = exp_diff_path.parts[-2]
            if comp_name not in self.comp_names:
                self.comp_names.append(comp_name)
            logging.info(f"Found {str(exp_diff_path)} of {comp_name}")

            if molecule_type in ["gene"]:
                self.gene_path_dic.setdefault(comp_name, exp_diff_path)
            elif molecule_type in ["isoform"]:
                self.tran_path_dic.setdefault(comp_name, exp_diff_path)

        self.comp_names.sort()

        return 1

    def stat_exp_diff(self, molecule_type):
        if molecule_type in ["gene"]:
            _diff_dic = self.gene_diff_dic
        elif molecule_type in ["isoform"]:
            _diff_dic = self.tran_diff_dic

        _stat_dic = dict()
        for _id, _sub_diff_dic in _diff_dic.items():
            for comp_name, info_dic in _sub_diff_dic.items():

                _stat_dic.setdefault(comp_name, {}).setdefault('DEG', 0)
                _stat_dic.setdefault(comp_name, {}).setdefault('DEG_Up', 0)
                _stat_dic.setdefault(comp_name, {}).setdefault('DEG_Down', 0)
                if info_dic['UpDown'] in ['up'] and info_dic['IsDEG'] in ['yes']:
                    _stat_dic[comp_name]['DEG'] += 1
                    _stat_dic[comp_name]['DEG_Up'] += 1
                elif info_dic['UpDown'] in ['down'] and info_dic['IsDEG'] in ['yes']:
                    _stat_dic[comp_name]['DEG'] += 1
                    _stat_dic[comp_name]['DEG_Down'] += 1

                _stat_dic.setdefault(comp_name, {}).setdefault('Delta_Pstv', 0)
                _stat_dic.setdefault(comp_name, {}).setdefault('Delta_Ngtv', 0)
                if info_dic['UpDown'] in ['up']:
                    _stat_dic[comp_name]['Delta_Pstv'] += 1
                elif info_dic['UpDown'] in ['down']:
                    _stat_dic[comp_name]['Delta_Ngtv'] += 1

        if molecule_type in ["gene"]:
            self.gene_stat_dic = _stat_dic
        elif molecule_type in ["isoform"]:
            self.tran_stat_dic = _stat_dic
        return 1


    def parse_exp_diff_path(self, molecule_type):
        if molecule_type in ["gene"]:
            _path_dic = self.gene_path_dic
        elif molecule_type in ["isoform"]:
            _path_dic = self.tran_path_dic

        _info_dic = dict()
        _diff_dic = dict()
        for comp_name, exp_diff_path in _path_dic.items():
            for line in exp_diff_path.open():
                items = line.rstrip('\n').split('\t')
                if items[0] in ["test_id"]:
                    idx_dic = dict()
                    for idx, item in enumerate(items):
                        idx_dic.setdefault(item, idx)
                    continue
                _id = items[idx_dic["test_id"]]
                _name = items[idx_dic["gene"]]
                _locus = items[idx_dic["locus"]]
                ctrl_id = items[idx_dic['sample_1']]
                case_id = items[idx_dic['sample_2']]
                status = items[idx_dic['status']]
                ctrl_value = items[idx_dic['value_1']]
                case_value = items[idx_dic['value_2']]
                log2fc = math.log(float(case_value)+1.0,2) - math.log(float(ctrl_value)+1.0,2)
                p = items[idx_dic['p_value']]
                q = items[idx_dic['q_value']]

                _info_dic.setdefault(_id, {}).setdefault("GeneName", _name)
                _info_dic.setdefault(_id, {}).setdefault("Locus", _locus)

                _diff_dic.setdefault(_id, {}).setdefault(comp_name, {})
                _diff_dic[_id][comp_name].setdefault('Status', status)
                _diff_dic[_id][comp_name].setdefault('CtrlName', ctrl_id)
                _diff_dic[_id][comp_name].setdefault('CaseName', case_id)
                _diff_dic[_id][comp_name].setdefault('CtrlFPKM', ctrl_value)
                _diff_dic[_id][comp_name].setdefault('CaseFPKM', case_value)
                _diff_dic[_id][comp_name].setdefault('Log2FC', log2fc)
                if log2fc > 0.0:
                    _diff_dic[_id][comp_name].setdefault('UpDown', 'up')
                elif log2fc < 0.0:
                    _diff_dic[_id][comp_name].setdefault('UpDown', 'down')
                else:
                    _diff_dic[_id][comp_name].setdefault('UpDown', 'flat')
                _diff_dic[_id][comp_name].setdefault('p', p)
                _diff_dic[_id][comp_name].setdefault('q', q)

                if self.p_criteria in [1.0] and self.q_criteria in [1.0]: # not use p, q
                    if abs(log2fc) >= self.log2fc_criteria:
                        _diff_dic[_id][comp_name].setdefault('IsDEG', 'yes')
                    else:
                        _diff_dic[_id][comp_name].setdefault('IsDEG', 'no')
                elif not self.p_criteria in [1.0] and self.q_criteria in [1.0]: # use p
                    if abs(log2fc) >= self.log2fc_criteria and float(p) < self.p_criteria:
                        _diff_dic[_id][comp_name].setdefault('IsDEG', 'yes')
                    else:
                        _diff_dic[_id][comp_name].setdefault('IsDEG', 'no')
                elif self.p_criteria in [1.0] and not self.q_criteria in [1.0]: # use q
                    if abs(log2fc) >= self.log2fc_criteria and float(q) < self.q_criteria:
                        _diff_dic[_id][comp_name].setdefault('IsDEG', 'yes')
                    else:
                        _diff_dic[_id][comp_name].setdefault('IsDEG', 'no')
                else: # usually... this way is worked.
                    if abs(log2fc) >= self.log2fc_criteria and float(p) < self.p_criteria and float(q) < self.q_criteria:
                        _diff_dic[_id][comp_name].setdefault('IsDEG', 'yes')
                    else:
                        _diff_dic[_id][comp_name].setdefault('IsDEG', 'no')

        if molecule_type in ["gene"]:
            self.gene_info_dic = _info_dic
            self.gene_diff_dic = _diff_dic
        elif molecule_type in ["isoform"]:
            self.tran_info_dic = _info_dic
            self.tran_diff_dic = _diff_dic

        return 1

    def anno_xxxx2desc(self, molecule_type, xxxx2desc_fp, colname, key):
        for line in xxxx2desc_fp.open():
            items = line.rstrip('\n').split('\t')
            if items[0] in ["ID"]:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            _id = items[idx_dic["ID"]]
            if molecule_type in ["gene"] and _id in self.gene_info_dic:
                self.gene_info_dic[_id].setdefault(key, items[idx_dic[colname]])
            if molecule_type in ["isoform"] and _id in self.tran_info_dic:
                self.tran_info_dic[_id].setdefault(key, items[idx_dic[colname]])

        if molecule_type in ["gene"] and key not in self.gene_info_col_s:
            self.gene_info_col_s.append(key)
        if molecule_type in ["isoform"] and key not in self.tran_info_col_s:
            self.tran_info_col_s.append(key)

        return 1

    def anno_xxxx2go(self, molecule_type, xxxx2go_fp):
        if molecule_type in ["gene"]:
            _info_dic = self.gene_info_dic
            _info_col_s = self.gene_info_col_s
        elif molecule_type in ["isoform"]:
            _info_dic = self.tran_info_dic
            _info_col_s = self.tran_info_col_s

        # add GO to headers
        for go_category in ["CC", "MF", "BP"]:
            if f"GO:{go_category}" not in _info_col_s:
                _info_col_s.append(f"GO:{go_category}")

        for line in xxxx2go_fp.open():
            items = line.rstrip('\n').split('\t')
            if items[0] in ["ID"]:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            _id = items[idx_dic["ID"]]
            go_id = items[idx_dic["ANNO_ID"]]
            _go_desc = items[idx_dic["ANNO_DESC"]]
            _go_category = _go_desc.split(':')[0]
            go_category = None
            if _go_category in ["C"]:
                go_category = "CC"
            elif _go_category in ["F"]:
                go_category = "MF"
            elif _go_category in ["P"]:
                go_category = "BP"
            go_term = _go_desc.split(':')[1]

            # init GO dic
            if _id in _info_dic:
                _info_dic[_id].setdefault("GO:CC", {})
                _info_dic[_id].setdefault("GO:MF", {})
                _info_dic[_id].setdefault("GO:BP", {})
                _info_dic[_id][f"GO:{go_category}"].setdefault(go_id, go_term)
            else:
                logging.warning(f"WARNING: {_id} is not in {molecule_type} info dic")


        if molecule_type in ["gene"]:
            self.gene_info_dic = _info_dic
            self.gene_info_col_s = _info_col_s
        elif molecule_type in ["isoform"]:
            self.tran_info_dic = _info_dic
            self.tran_info_col_s = _info_col_s

        return 1

    def write_deg_tsv(self, molecule_type, _outfn):
        if molecule_type in ["gene"]:
            _info_dic = self.gene_info_dic
            _diff_dic = self.gene_diff_dic
            _info_col_s = self.gene_info_col_s
        elif molecule_type in ["isoform"]:
            _info_dic = self.tran_info_dic
            _diff_dic = self.tran_diff_dic
            _info_col_s = self.tran_info_col_s

        headers = ["GeneID"]
        for comp_name in self.comp_names:
            for _diff_col in self._diff_col_s:
                headers.append('{0}:{1}'.format(comp_name, _diff_col))
        headers.extend(_info_col_s)

        outfp = Path(f"{self.outprefix}.{_outfn}")
        logging.info(f"writing DEG table... {str(outfp)}")
        outfh = outfp.open("w")
        outfh.write("{0}\n".format("\t".join([str(x) for x in headers])))
        for _id, _sub_dic in _diff_dic.items():
            items = [_id]
            for comp_name in self.comp_names:
                for _diff_col in self._diff_col_s:
                    if _diff_col in ["CtrlFPKM", "CaseFPKM", "Log2FC"]:
                        items.append(round(float(_sub_dic[comp_name][_diff_col]), 4))
                    else:
                        items.append(_sub_dic[comp_name][_diff_col])
            for _info_col in _info_col_s:
                if _info_col in _info_dic[_id]:
                    items.append(_info_dic[_id][_info_col])
                else:
                    items.append("None")
            outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
        outfh.close()
        return 1

    def write_deg_stat_tsv(self, molecule_type, _outfn):
        if molecule_type in ["gene"]:
            _stat_dic = self.gene_stat_dic
        elif molecule_type in ["isoform"]:
            _stat_dic = self.tran_stat_dic

        headers = ["ComparisonID"]
        headers.append("DEG")
        headers.append("DEG_Up")
        headers.append("DEG_Down")
        headers.append("Delta_Pstv")
        headers.append("Delta_Ngtv")

        outfp = Path(f"{self.outprefix}.{_outfn}")
        logging.info(f"writing DEG Stat table... {str(outfp)}")
        outfh = outfp.open("w")
        outfh.write("{0}\n".format("\t".join([str(x) for x in headers])))
        for comp_name in self.comp_names:
            info_dic = _stat_dic[comp_name]
            items = [comp_name]
            items.append(info_dic["DEG"])
            items.append(info_dic["DEG_Up"])
            items.append(info_dic["DEG_Down"])
            items.append(info_dic["Delta_Pstv"])
            items.append(info_dic["Delta_Ngtv"])
            outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
        outfh.close()

        return 1


class STAR:
    def __init__(self, wkdir, outprefix, target_dir):
        self.wkdir = Path(wkdir)
        self.input_read_type = None

        # parse Log.final.out files
        self.align_summary_path_dic = dict()
        self.align_summary_info_dic = dict()
        self.find_align_summary_path(target_dir)
        self.parse_align_summary_path()
        self.write_align_summary_tsv(Path(f"{outprefix}.xls"))

    def write_align_summary_tsv(self, file_path):
        items = ["Samples"]
        items.append("InputReads")
        items.append("MappedReads")
        items.append("MappedReadsUniquely")
        items.append("MappedReadsMultiple")
        ofh = file_path.open("w")
        ofh.write("{0}\n".format("\t".join(items)))
        for sample_name, _dic in sorted(self.align_summary_info_dic.items()):
            items = [sample_name]
            items.append(_dic["input reads"])
            items.append(_dic["mapped reads"])
            items.append(_dic["uniquely mapped reads"])
            items.append(_dic["multiple mapped reads"])
            ofh.write("{0}\n".format("\t".join([str(x) for x in items])))
        ofh.close()
        return 1


    def find_align_summary_path(self, _dir):
        for align_summary_path in Path(_dir).rglob("*/Log.final.out"):
            sample_name = align_summary_path.parts[-2]
            logging.info(f"Found {str(align_summary_path)} of {sample_name}")
            self.align_summary_path_dic.setdefault(sample_name, align_summary_path)
        return 1

    def parse_align_summary_path(self):
        for sample_name, align_summary_path in self.align_summary_path_dic.items():
            _dic = self.star_logFinal_parser(align_summary_path)
            self.align_summary_info_dic.setdefault(sample_name, _dic)

    def star_logFinal_parser(self, file_path):
        in_read = 0
        mapped_reads = 0
        mapped_reads_uniquely = 0
        mapped_reads_multiple = 0
        unmapped_reads = 0
        for line in file_path.open():
            if not line.strip():
                continue
            items = line.strip().split('\t')
            if items[0].startswith('Number of input reads'):
                in_read = int(items[1])*2
            elif items[0].startswith('Uniquely mapped reads number'):
                mapped_reads += int(items[1])*2
                mapped_reads_uniquely = int(items[1])*2
            elif items[0].startswith('Number of reads mapped to multiple loci'):
                mapped_reads += int(items[1])*2
                mapped_reads_multiple += int(items[1])*2
            elif items[0].startswith('Number of reads mapped to too many loci'):
                mapped_reads += int(items[1])*2
                mapped_reads_multiple += int(items[1])*2
            elif items[0].startswith('Number of reads unmapped: too many mismatches'):
                unmapped_reads += int(items[1])*2
            elif items[0].startswith('Number of reads unmapped: too short'):
                unmapped_reads += int(items[1])*2
            elif items[0].startswith('Number of reads unmapped: other'):
                unmapped_reads += int(items[1])*2

        if not mapped_reads in [mapped_reads_uniquely+mapped_reads_multiple]:
            print('ERROR: check the star_logFinal_parser')
            print('if not mapped_reads in [mapped_reads_uniquely+mapped_reads_multiple]:')
            print('in_read', in_read)
            print('mapped_reads', mapped_reads)
            print('mapped_reads_uniquely', mapped_reads_uniquely)
            print('mapped_reads_multiple', mapped_reads_multiple)
            print('unmapped_reads', unmapped_reads)
            sys.exit()
        if not in_read in [mapped_reads+unmapped_reads]:
            print('ERROR: check the star_logFinal_parser')
            print('if not in_read in [mapped_reads+unmapped_reads]:')
            print('in_read', in_read)
            print('mapped_reads', mapped_reads)
            print('mapped_reads_uniquely', mapped_reads_uniquely)
            print('mapped_reads_multiple', mapped_reads_multiple)
            print('unmapped_reads', unmapped_reads)
            sys.exit()

        _dic = dict()
        _dic.setdefault('input reads', in_read)
        _dic.setdefault('mapped reads', mapped_reads)
        _dic.setdefault('unmapped reads', unmapped_reads)
        _dic.setdefault('uniquely mapped reads', mapped_reads_uniquely)
        _dic.setdefault('multiple mapped reads', mapped_reads_multiple)
        return _dic


class Rsem:
    def __init__(self, wkdir, outprefix, gene2desc, tran2desc, sample_id_s):
        self.rsem_dic = dict()
        self.gene_info_dic = dict()
        self.tran_info_dic = dict()

        self.wkdir = Path(wkdir)
        self.outprefix = outprefix
        self.sample_id_s = sample_id_s

        self.gene2desc_fp = Path(gene2desc)
        self.tran2desc_fp = Path(tran2desc)

    def anno_xxxx2desc(self, molecule_type, xxxx2desc_fp, colname, key):
        if molecule_type in ["genes"]:
            _info_dic = self.gene_info_dic
        elif molecule_type in ["isoforms"]:
            _info_dic = self.tran_info_dic

        for line in xxxx2desc_fp.open():
            items = line.rstrip('\n').split('\t')
            if items[0] in ["ID"]:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            _id = items[idx_dic["ID"]]
            _info_dic.setdefault(_id, {}).setdefault(key, items[idx_dic[colname]])

        if molecule_type in ["genes"]:
            self.gene_info_dic = _info_dic
        elif molecule_type in ["isoforms"]:
            self.tran_info_dic = _info_dic
        return 1

    def dic_to_stat(self, molecule_type):
        if molecule_type in ["genes"]:
            _info_dic = self.gene_info_dic
        elif molecule_type in ["isoforms"]:
            _info_dic = self.tran_info_dic

        _stat_col_s = ["fpkm>=0.3", "fpkm>=1", "fpkm>=10", "fpkm>=100", "fpkm>=1000"]

        stat_dic = dict()
        for sample_id in self.sample_id_s:
            for _stat_col in _stat_col_s:
                stat_dic.setdefault(sample_id, {}).setdefault(_stat_col, [])

        for _id, a_dic in self.rsem_dic.items():
            biotype = _info_dic[_id]['biotype']
            if biotype not in ['protein_coding']:
                continue
            for sample_id in self.sample_id_s:
                fpkm = float(a_dic[sample_id]['fpkm'])
                if fpkm >= 1000:
                    stat_dic[sample_id]["fpkm>=1000"].append(_id)
                if fpkm >= 100:
                    stat_dic[sample_id]["fpkm>=100"].append(_id)
                if fpkm >= 10:
                    stat_dic[sample_id]["fpkm>=10"].append(_id)
                if fpkm >= 1:
                    stat_dic[sample_id]["fpkm>=1"].append(_id)
                if fpkm >= 0.3:
                    stat_dic[sample_id]["fpkm>=0.3"].append(_id)

        outfn = f"{self.outprefix}.{molecule_type}.stats.exp.xls"
        outfh = Path(outfn).open("w")
        headers = ["SampleNames"]
        headers.extend(_stat_col_s)
        outfh.write("{0}\n".format("\t".join(headers)))
        for sample_id in self.sample_id_s:
            items = [sample_id]
            for _stat_col in _stat_col_s:
                items.append(len(stat_dic[sample_id][_stat_col]))
            outfh.write("{0}\n".format("\t".join([str(x) for x in items])))
        outfh.close()

    def dic_to_table(self, _type, data, threshold=0.0):
        outfn = f"{self.outprefix}.{_type}.{data}.xls"
        outfh = Path(outfn).open("w")

        cols = ['ID']
        cols.extend(self.sample_id_s)
        cols.extend(["Locus","GeneName","GeneSymbol","Description","Biotype"])
        outfh.write("{0}\n".format('\t'.join(cols)))

        for _id, s_dic in sorted(self.rsem_dic.items()):
            items = [_id]
            values = []
            for s_id in self.sample_id_s:
                if data in ["fpkm.all", "fpkm.exp"]:
                    value = s_dic[s_id]["fpkm"]
                else:
                    value = s_dic[s_id][data]
                items.append(value)
                values.append(value)

            if _type in ["genes"] and _id in self.gene_info_dic:
                items.append("")
                items.append(self.gene_info_dic[_id]['symbol'])
                items.append(self.gene_info_dic[_id]['symbol'])
                items.append(self.gene_info_dic[_id]['description'])
                items.append(self.gene_info_dic[_id]['biotype'])
            elif _type in ["isoforms"] and _id in self.tran_info_dic:
                items.append("")
                items.append(self.tran_info_dic[_id]['symbol'])
                items.append(self.tran_info_dic[_id]['symbol'])
                items.append(self.tran_info_dic[_id]['description'])
                items.append(self.tran_info_dic[_id]['biotype'])
            else:
                items.append("-")
                items.append("-")
                items.append("-")

            is_write = 1
            values = [float(x) for x in values]
            if max(values) < threshold:
                is_write = 0

            if is_write:
                outfh.write('{0}\n'.format('\t'.join(items)))

        outfh.close()

    def from_fn_to_dic(self, _type):
        for sample_id in self.sample_id_s:
            _fp = self.wkdir / "anl" / "rsem" / sample_id / "rsem.{0}.results".format(_type)
            if _type in ['genes']:
                self.parse_gene_fn(sample_id, _fp)
            elif _type in ['isoforms']:
                self.parse_isoforms_fn(sample_id, _fp)

    def parse_gene_fn(self, sid, gene_fp):
        for line in gene_fp.open():
            items = line.rstrip('\n').split('\t')
            if items[0] in ['gene_id']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            _id = items[idx_dic['gene_id']]
            length = str(int(round(float(items[idx_dic['length']]),0)))
            readcount = str(int(round(float(items[idx_dic['expected_count']]),0)))
            tpm = items[idx_dic['TPM']]
            fpkm = items[idx_dic['FPKM']]

            self.rsem_dic.setdefault(_id,{}).setdefault(sid,{})
            self.rsem_dic[_id][sid].setdefault('length', length)
            self.rsem_dic[_id][sid].setdefault('readcount', readcount)
            self.rsem_dic[_id][sid].setdefault('tpm', tpm)
            self.rsem_dic[_id][sid].setdefault('fpkm', fpkm)

    def parse_isoforms_fn(self, sid, isoform_fp):
        for line in isoform_fp.open():
            items = line.rstrip('\n').split('\t')
            if items[0] in ['transcript_id']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            _id = items[idx_dic['transcript_id']]
            length = str(int(round(float(items[idx_dic['length']]),0)))
            readcount = str(int(round(float(items[idx_dic['expected_count']]),0)))
            tpm = items[idx_dic['TPM']]
            fpkm = items[idx_dic['FPKM']]

            self.rsem_dic.setdefault(_id,{}).setdefault(sid,{})
            self.rsem_dic[_id][sid].setdefault('length', length)
            self.rsem_dic[_id][sid].setdefault('readcount', readcount)
            self.rsem_dic[_id][sid].setdefault('tpm', tpm)
            self.rsem_dic[_id][sid].setdefault('fpkm', fpkm)




def main(args):
    logging.info(args)

    if args.biotool in ['fastp']:

        fastp = Fastp(args.wkdir, args.outprefix)
        fastp.find_json_path(args.target_dir)
        fastp.write_summary_tsv(args.target_tag)

    elif args.biotool in ['tophat']:

        tophat = Tophat(args.wkdir, args.outprefix, args.target_dir)

    elif args.biotool in ['star']:

        star= STAR(args.wkdir, args.outprefix, args.target_dir)

    elif args.biotool in ['cufflinks']:

        cufflinks = Cufflinks(args.wkdir, args.outprefix, args.gene2desc, args.tran2desc, args.gene2go, args.tran2go)

        cufflinks.find_fpkm_tracking_path("genes", "anl/cufflinks")
        if cufflinks.gene2desc_fp.is_file():
            cufflinks.anno_xxxx2desc("genes", cufflinks.gene2desc_fp, "BIOTYPE", "biotype")
            cufflinks.anno_xxxx2desc("genes", cufflinks.gene2desc_fp, "SYMBOL", "symbol")
            cufflinks.anno_xxxx2desc("genes", cufflinks.gene2desc_fp, "DESCRIPTION", "description")
        cufflinks.parse_fpkm_tracking_path("genes")
        cufflinks.write_fpkm_tsv("genes", 0.0, cufflinks.wkdir / f"{args.outprefix}.genes.fpkm.all.xls")
        cufflinks.write_fpkm_tsv("genes", 0.3, cufflinks.wkdir / f"{args.outprefix}.genes.fpkm.exp.xls")
        cufflinks.write_stat_exp("genes", cufflinks.wkdir / f"{args.outprefix}.genes.stats.exp.xls")
        cufflinks.write_stat_coexp("genes", cufflinks.wkdir / f"{args.outprefix}.genes.stats.coexp.xls")

        cufflinks.find_fpkm_tracking_path("isoforms", "anl/cufflinks")
        if cufflinks.tran2desc_fp.is_file():
            cufflinks.anno_xxxx2desc("isoforms", cufflinks.tran2desc_fp, "BIOTYPE", "biotype")
            cufflinks.anno_xxxx2desc("isoforms", cufflinks.tran2desc_fp, "SYMBOL", "symbol")
            cufflinks.anno_xxxx2desc("isoforms", cufflinks.tran2desc_fp, "DESCRIPTION", "description")
        cufflinks.parse_fpkm_tracking_path("isoforms")
        cufflinks.write_fpkm_tsv("isoforms", 0.0, cufflinks.wkdir / f"{args.outprefix}.isoforms.fpkm.all.xls")
        cufflinks.write_fpkm_tsv("isoforms", 0.3, cufflinks.wkdir / f"{args.outprefix}.isoforms.fpkm.exp.xls")
        cufflinks.write_stat_exp("isoforms", cufflinks.wkdir / f"{args.outprefix}.isoforms.stats.exp.xls")
        cufflinks.write_stat_coexp("isoforms", cufflinks.wkdir / f"{args.outprefix}.isoforms.stats.coexp.xls")


    elif args.biotool in ['cuffnorm']:

        cuffnorm = Cuffnorm(args.wkdir, args.outprefix, args.gene2desc, args.tran2desc, args.gene2go, args.tran2go)

        cuffnorm.find_fpkm_table_path("genes", "anl/cuffnorm")
        if cuffnorm.gene2desc_fp.is_file():
            cuffnorm.anno_xxxx2desc("genes", cuffnorm.gene2desc_fp, "BIOTYPE", "biotype")
            cuffnorm.anno_xxxx2desc("genes", cuffnorm.gene2desc_fp, "SYMBOL", "symbol")
            cuffnorm.anno_xxxx2desc("genes", cuffnorm.gene2desc_fp, "DESCRIPTION", "description")
        if cuffnorm.gene2go_fp.is_file():
            cuffnorm.anno_xxxx2go("genes", cuffnorm.gene2go_fp)
        cuffnorm.parse_fpkm_table_path("genes")
        cuffnorm.write_fpkm_tsv("genes", 0.0, cuffnorm.wkdir / f"{args.outprefix}.genes.fpkm.all.xls")
        cuffnorm.write_fpkm_tsv("genes", 0.3, cuffnorm.wkdir / f"{args.outprefix}.genes.fpkm.exp.xls")
        cuffnorm.write_stat_exp("genes", cuffnorm.wkdir / f"{args.outprefix}.genes.stats.exp.xls")
        cuffnorm.write_stat_coexp("genes", cuffnorm.wkdir / f"{args.outprefix}.genes.stats.coexp.xls")

        cuffnorm = Cuffnorm(args.wkdir, args.outprefix, args.gene2desc, args.tran2desc, args.gene2go, args.tran2go)

        cuffnorm.find_fpkm_table_path("isoforms", "anl/cuffnorm")
        if cuffnorm.tran2desc_fp.is_file():
            cuffnorm.anno_xxxx2desc("isoforms", cuffnorm.tran2desc_fp, "BIOTYPE", "biotype")
            cuffnorm.anno_xxxx2desc("isoforms", cuffnorm.tran2desc_fp, "SYMBOL", "symbol")
            cuffnorm.anno_xxxx2desc("isoforms", cuffnorm.tran2desc_fp, "DESCRIPTION", "description")
        cuffnorm.parse_fpkm_table_path("isoforms")
        cuffnorm.write_fpkm_tsv("isoforms", 0.0, cuffnorm.wkdir / f"{args.outprefix}.isoforms.fpkm.all.xls")
        cuffnorm.write_fpkm_tsv("isoforms", 0.3, cuffnorm.wkdir / f"{args.outprefix}.isoforms.fpkm.exp.xls")
        cuffnorm.write_stat_exp("isoforms", cuffnorm.wkdir / f"{args.outprefix}.isoforms.stats.exp.xls")
        cuffnorm.write_stat_coexp("isoforms", cuffnorm.wkdir / f"{args.outprefix}.isoforms.stats.coexp.xls")

    elif args.biotool in ['rsem']:
        rsem = Rsem(args.wkdir, args.outprefix, args.gene2desc, args.tran2desc, args.sample_id_s)
        if rsem.gene2desc_fp.is_file():
            rsem.anno_xxxx2desc("genes", rsem.gene2desc_fp, "BIOTYPE", "biotype")
            rsem.anno_xxxx2desc("genes", rsem.gene2desc_fp, "SYMBOL", "symbol")
            rsem.anno_xxxx2desc("genes", rsem.gene2desc_fp, "DESCRIPTION", "description")
        if rsem.tran2desc_fp.is_file():
            rsem.anno_xxxx2desc("isoforms", rsem.tran2desc_fp, "BIOTYPE", "biotype")
            rsem.anno_xxxx2desc("isoforms", rsem.tran2desc_fp, "SYMBOL", "symbol")
            rsem.anno_xxxx2desc("isoforms", rsem.tran2desc_fp, "DESCRIPTION", "description")
        rsem.rsem_dic = dict()
        rsem.from_fn_to_dic('genes')
        rsem.dic_to_table('genes', 'length')
        rsem.dic_to_table('genes', 'readcount')
        rsem.dic_to_table('genes', 'tpm')
        rsem.dic_to_table('genes', 'fpkm.all')
        rsem.dic_to_table('genes', 'fpkm.exp', 0.3)
        rsem.dic_to_stat('genes')
        rsem.rsem_dic = dict()
        rsem.from_fn_to_dic('isoforms')
        rsem.dic_to_table('isoforms', 'length')
        rsem.dic_to_table('isoforms', 'readcount')
        rsem.dic_to_table('isoforms', 'tpm')
        rsem.dic_to_table('isoforms', 'fpkm.all')
        rsem.dic_to_table('isoforms', 'fpkm.exp', 0.3)
        rsem.dic_to_stat('isoforms')

    elif args.biotool in ['cuffdiff']:

        cuffdiff = Cuffdiff(args.wkdir, args.outprefix, args.log2fc_criteria, args.p_criteria, args.q_criteria)

        if Path(args.gene2desc).is_file():
            cuffdiff.anno_xxxx2desc("gene", Path(args.gene2desc), "SYMBOL", "GeneSymbol")
            cuffdiff.anno_xxxx2desc("gene", Path(args.gene2desc), "DESCRIPTION", "Description")
            cuffdiff.anno_xxxx2desc("gene", Path(args.gene2desc), "BIOTYPE", "Biotype")

        if Path(args.tran2desc).is_file():
            cuffdiff.anno_xxxx2desc("isoform", Path(args.tran2desc), "SYMBOL", "TranscriptSymbol")
            cuffdiff.anno_xxxx2desc("isoform", Path(args.tran2desc), "DESCRIPTION", "Description")
            cuffdiff.anno_xxxx2desc("isoform", Path(args.tran2desc), "BIOTYPE", "Biotype")

        if Path(args.gene2go).is_file():
            cuffdiff.anno_xxxx2go("gene", Path(args.gene2go))
        if Path(args.tran2go).is_file():
            #cuffdiff.anno_xxxx2go("isoform", Path(args.tran2go))
            pass # file size is too large

        cuffdiff.write_deg_tsv("gene", "genes.xls")
        cuffdiff.write_deg_tsv("isoform", "isoforms.xls")
        cuffdiff.write_deg_stat_tsv("gene", "genes.stats.xls")
        cuffdiff.write_deg_stat_tsv("isoform", "isoforms.stats.xls")

    elif args.biotool in ['get_anno_uniprot']:
        get_anno_uniprot = GetAnnoUniprot(args.wkdir, args.gtf, args.tsv, args.dat)
    elif args.biotool in ['get_mapg2t_gtf']:
        get_mapg2t_gtf = GetMapG2TGTF(args.wkdir, args.gtf)
    elif args.biotool in ['get_desc']:
        get_desc = GetDesc(args.wkdir, args.gtf, args.mysqlgene, args.mysqltran)
    elif args.biotool in ['cleanup_fasta']:
        cleanup_fasta = CleanupFasta(args.wkdir, args.target, args.gtf, args.infasta, args.outfasta)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description="MultiParser")
    parser.add_argument("--wkdir", default="./")
    subparsers = parser.add_subparsers(required=True, dest="biotool")

    subparser = subparsers.add_parser("fastp")
    subparser.add_argument("--outprefix", default="MultiParser.Fastp.Raw") # MultiParser.Fastp.Clean if cleanfastq
    subparser.add_argument("--target-dir", default="readyfastq") # cleanfastq
    subparser.add_argument("--target-tag", default="after_filtering", choices=("after_filtering", "before_filtering"))

    subparser = subparsers.add_parser("tophat")
    subparser.add_argument("--outprefix", default="MultiParser.Tophat")
    subparser.add_argument("--target-dir", default="anl/tophat")

    subparser = subparsers.add_parser("star")
    subparser.add_argument("--outprefix", default="MultiParser.STAR")
    subparser.add_argument("--target-dir", default="anl/star")

    subparser = subparsers.add_parser("cufflinks")
    subparser.add_argument("--outprefix", default="MultiParser.Cufflinks")
    subparser.add_argument("--gene2desc", default="anno.gene2desc.tsv")
    subparser.add_argument("--tran2desc", default="anno.transcript2desc.tsv")

    subparser = subparsers.add_parser("cuffnorm")
    subparser.add_argument("--outprefix", default="MultiParser.Cuffnorm")
    subparser.add_argument("--gene2desc", default="anno.gene2desc.tsv")
    subparser.add_argument("--tran2desc", default="anno.transcript2desc.tsv")
    subparser.add_argument("--gene2go", default="anno.gene2go.tsv")
    subparser.add_argument("--tran2go", default="anno.transcript2go.tsv")

    subparser = subparsers.add_parser("rsem")
    subparser.add_argument("--outprefix", default="MultiParser.Rsem")
    subparser.add_argument("--gene2desc", default="anno.gene2desc.tsv")
    subparser.add_argument("--tran2desc", default="anno.transcript2desc.tsv")
    subparser.add_argument("--sample-id-s", nargs="+", default=[])

    subparser = subparsers.add_parser("cuffdiff")
    subparser.add_argument("--outprefix", default="MultiParser.Cuffdiff")
    subparser.add_argument("--gene2desc", default="anno.gene2desc.tsv")
    subparser.add_argument("--tran2desc", default="anno.transcript2desc.tsv")
    subparser.add_argument("--gene2go", default="anno.gene2go.tsv")
    subparser.add_argument("--tran2go", default="anno.transcript2go.tsv")
    subparser.add_argument("--log2fc-criteria", default=0.584, type=float) # fc 1.5
    subparser.add_argument("--p-criteria", default=0.05, type=float)
    subparser.add_argument("--q-criteria", default=1.0, type=float)

    subparser = subparsers.add_parser('get_anno_uniprot')
    subparser.add_argument('--gtf')
    subparser.add_argument('--tsv')
    subparser.add_argument('--dat')

    subparser = subparsers.add_parser('get_mapg2t_gtf')
    subparser.add_argument('--gtf')

    subparser = subparsers.add_parser('get_desc')
    subparser.add_argument('--gtf')
    subparser.add_argument('--mysqlgene')
    subparser.add_argument('--mysqltran')

    subparser = subparsers.add_parser('cleanup_fasta')
    subparser.add_argument('--target')
    subparser.add_argument('--gtf')
    subparser.add_argument('--infasta')
    subparser.add_argument('--outfasta')

    args = parser.parse_args()
    main(args)






