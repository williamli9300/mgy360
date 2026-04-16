from pathlib import Path

output_lines = []

path = Path("./")

for file in path.iterdir():
    if ".vcf" in file.name:
        with open(file.name, "r") as f:
            lines = f.readlines()
        id_ = ""
        for line in lines:
            if len(line) > 1:
                if line[1] != "#":
                    current_line = line.split("\t")
                    if current_line[0] == "#CHROM":
                        id_ = current_line[9].replace("\n", "")
                    elif "INDEL" not in current_line[7]:
                        current_line[-1] = current_line[-1].replace("\n", "")
                        current_line.append(id_)
                        output_lines.append("\t".join(current_line))

with open("./master.vcf", "w") as f:
    f.writelines(f"{line}\n" for line in output_lines)
