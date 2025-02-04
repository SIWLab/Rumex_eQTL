f_in = open("noMatPAR_txanno.siteTypes", "r")
out_file1 = "nonsyn.sites"
out_file2 = "syn.sites"

f_out1 = open(out_file1, "w")
f_out2 = open(out_file2, "w")

content = f_in.readlines()
for line in content:
    site_line = line.split()
    degeneracy = site_line[4]
    if degeneracy == "0":
        print(line.strip(), file=f_out1)
    elif degeneracy == "4":
        print(line.strip(), file=f_out2)

f_in.close()
f_out1.close()
f_out2.close