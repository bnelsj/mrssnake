import os
import argparse

import pandas as pd

def get_samples_for_mapping(watch_dir, manifest_size):
    mapped = []
    for fn in os.listdir(watch_dir + "/" + "mapping_finished"):
        sn = fn.replace(".txt", "")
        if sn not in mapped:
            mapped.append(sn)

    unmapped = []
    for fn in os.listdir(watch_dir + "/" + "download_finished"):
        sn = fn.replace(".txt", "")
        if sn not in mapped and sn not in unmapped:
            unmapped.append(sn)
    
    to_map = []
    for fn in unmapped:
        file = "%s/%s/%s.txt" % (watch_dir, "currently_mapping", fn)
        if not os.path.exists(file):
            open(file, "a").close()
            to_map.append(fn)
            if len(to_map) == manifest_size:
                break

    return to_map

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("watch_dir", help="Directory with 'download_finished', 'mapping_finished' \
                        and 'bam' subdirectories")
    parser.add_argument("manifest", help="Path to output manifest file")
    parser.add_argument("manifest_size", type=int, default=10, help="Number of samples to include in output manifest")

    args = parser.parse_args()

    samples = []

    to_map = get_samples_for_mapping(args.watch_dir, args.manifest_size)

    dat = pd.DataFrame(data={"sn": to_map, "source": "local", "bam": "", "index": ""})

    dat["bam"] = dat.sn.map(lambda x: "%s/bam/%s.bam" % (args.watch_dir, x))
    dat["index"] = dat.sn.map(lambda x: "%s/bam/%s.bai" % (args.watch_dir, x))

    dat.to_csv(args.manifest, sep="\t", index=False)
