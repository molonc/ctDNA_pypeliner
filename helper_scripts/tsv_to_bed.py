import csv

def main():
    tsv_file = open('beds/CG001.v3.4_Hotspot_Manifest_Panel3.4.5_20170921.tsv', 'r')
    bed_file = open('beds/CG001.v3.4', 'w+')
    tsv_reader = csv.DictReader(v_four, delimiter="\t")
    bed_writer = csv.writer(v_four_bed, delimiter="\t")
    for row in tsv_reader:
        bed_writer.writerow([row['Chr'], row['Start'], row['End']])
        bed_file.flush()



if __name__ == '__main__':
    main()