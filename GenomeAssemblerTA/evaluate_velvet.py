import main
import os
from contigs3 import calc_n50, longest_contig

def run():
    velvet_path = r"../data/results/"
    summary_path = "../data/results/velvet_summary.txt"
    main.write_file(summary_path, ("PATH", "N50", "LONGEST", "# of contigs"))

    for f in os.listdir(velvet_path):
        if f != "old" and os.path.isdir(os.path.join(velvet_path,f)):
            new_path = os.path.join(os.path.join(velvet_path, f), "contigs.fa")
            res = main.open_file(new_path)
            res = "".join(["," if ">NODE" in x else x for x in res]).split(",") # combine nodes etc.
            res = [x for x in res if x!='']
            print(res)
            n50 = calc_n50(res)
            lng = longest_contig(res)
            summary = (f, n50, lng, len(res))
            main.write_file(summary_path,summary, mode="a+",newline=True)



if __name__ == '__main__':
    run()