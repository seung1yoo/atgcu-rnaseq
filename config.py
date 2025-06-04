

import yaml
from pathlib import Path
import pandas as pd


class ConfigMaker:
    def __init__(self, config_dir):
        Path(config_dir).mkdir(exist_ok=True)
        self.sample_conf_dic = {"samples": {}}
        self.compare_conf_dic = {"compares": {}}

    def update_sample_conf_dic(self, _path, _tag, _fastq_key):
        sample_name = _path.name.split(_tag)[0]
        if sample_name in self.sample_conf_dic["samples"]:
            self.sample_conf_dic["samples"][sample_name][_fastq_key] = str(_path)
        else:
            print(f"WARNNING : {sample_name} {str(_path)} is not in meta/sample_meta.csv. plz double check!")

        return 1

    def load_sample_meta(self, _path):
        for line in _path.open():
            items = line.rstrip().split(',')
            if items[0] in ['sample']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            sample_name = items[idx_dic["sample"]]
            group = items[idx_dic["group"]]

            self.sample_conf_dic["samples"].setdefault(sample_name, {})
            for _fastq_key in ["fastq_1", "fastq_2"]:
                self.sample_conf_dic["samples"][sample_name].setdefault(_fastq_key, Path())
            self.sample_conf_dic["samples"][sample_name].setdefault("cuff_librarytype", "fr-unstranded")
            self.sample_conf_dic["samples"][sample_name].setdefault("group", group)

    def grep_fastq_gz(self, fastq_dir):
        for _path in Path(fastq_dir).glob("*.fastq.gz"):
            if _path.name.endswith("_R1.fastq.gz"):
                self.update_sample_conf_dic(_path, "_R1.fastq.gz", "fastq_1")
            elif _path.name.endswith("_1.fastq.gz"):
                self.update_sample_conf_dic(_path, "_1.fastq.gz", "fastq_1")
            elif _path.name.endswith("_R2.fastq.gz"):
                self.update_sample_conf_dic(_path, "_R2.fastq.gz", "fastq_2")
            elif _path.name.endswith("_2.fastq.gz"):
                self.update_sample_conf_dic(_path, "_2.fastq.gz", "fastq_2")
            else:
                continue
        return 1

    def write_sample_conf(self, outfn):
        with open(outfn, "w") as outfh:
            yaml.dump(self.sample_conf_dic, outfh)
        return 1

    def grep_group_samples(self, group_id):
        samples = list()
        for sample_name, info_dic in self.sample_conf_dic["samples"].items():
            if info_dic["group"] in [group_id]:
                samples.append(sample_name)
        return samples


    def prepare_compare(self, _path):
        for line in _path.open():
            items = line.rstrip().split(',')
            if items[0] in ['compid']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            compare_id = items[idx_dic['compid']]
            control = items[idx_dic['control']]
            case = items[idx_dic['case']]

            control_samples = self.grep_group_samples(control)
            case_samples = self.grep_group_samples(case)

            self.compare_conf_dic["compares"].setdefault(compare_id, {})
            self.compare_conf_dic["compares"][compare_id].setdefault("control", {}).setdefault("group_name", control)
            self.compare_conf_dic["compares"][compare_id].setdefault("control", {}).setdefault("samples", control_samples)
            self.compare_conf_dic["compares"][compare_id].setdefault("treatment", {}).setdefault("group_name", case)
            self.compare_conf_dic["compares"][compare_id].setdefault("treatment", {}).setdefault("samples", case_samples)

    def write_compare_conf(self, outfn):
        with open(outfn, "w") as outfh:
            yaml.dump(self.compare_conf_dic, outfh)
        return 1


def main(args):


    configmaker = ConfigMaker(args.config_dir)

    if Path(args.sample_meta).exists():
        configmaker.load_sample_meta(Path(args.sample_meta))
        configmaker.grep_fastq_gz(args.fastq_dir)
    configmaker.write_sample_conf(f"{args.config_dir}/sample.yml")

    if Path(args.compare_meta).exists():
        configmaker.prepare_compare(Path(args.compare_meta))
    configmaker.write_compare_conf(f"{args.config_dir}/compare.yml")




if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq-dir", default="readyfastq")
    parser.add_argument("--config-dir", default="config")
    parser.add_argument("--sample-meta", default="meta/sample_meta.csv")
    parser.add_argument("--compare-meta", default="meta/compare_meta.csv")
    args = parser.parse_args()
    main(args)
