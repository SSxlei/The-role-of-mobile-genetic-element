import os

try:
    os.mkdir("ACT")
except FileExistsError:
    print("Directory 'ACT' already exists.")

listfile = os.listdir(".")
for line in listfile:  
    if line.endswith(".fna"):
        file_name = line[:-4]  # Remove the ".fna" suffix
        oldfile = open(line, "r")
        newfile = open("ACT/" + file_name + ".fna", "w")
        temp_i = 0
        for cont in oldfile:
            if cont.startswith('>'):
                ids = '>' + file_name
                temp_i += 1
                newfile.write(ids + "\n")
            else:
                newfile.write(cont)
        newfile.close()
        oldfile.close()
