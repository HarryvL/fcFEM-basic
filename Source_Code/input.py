import sys

nodes=[]
elements=[]
elset=[]
nset=[]
sload=[]
bcon=[]
matset=[]

def data_lines(index):
    if index +1 < len(file_clean):
        while not file_clean[index+1][0].startswith('*'):
            index+=1
            print("data: ", file_clean[index])
    return index

def key_words(index):
    key_list = []
    for keyw in file_clean[index][1:]:
        key_list.append(keyw)
    return key_list

with open('../FEMMeshGmsh.inp', encoding='utf8') as f:
    file_raw = f.readlines()

file_clean = []

for line in file_raw:
    upper = line.upper()
    if upper.startswith("**") or upper.strip() == "":
        pass
    else:
        list = (upper.strip()).split(',')
        file_clean.append(list)

# for line in file_clean:
#     print(line)
# print("x0x0x0x0x0x0x")
# sys.exit()
for index in range(len(file_clean)):
    first = file_clean[index][0]
    if first.startswith("*"):
        print(file_clean[index])
        if first == "*HEADING":
            index = data_lines(index)
        elif first == "*RESTART":
            index = data_lines(index)
        elif first == "*NODE":
            key_list = key_words(index)
            print(key_list)
            print(file_clean[index+1])
            while not file_clean[index + 1][0].startswith('*'):
                index += 1
                a = file_clean[index]
                nodes.append([int(a[0]), float(a[1]), float(a[2]), float(a[3])])
        elif first == "*NODE FILE":
            key_list = key_words(index)
            index = data_lines(index)
        elif first == "*ELEMENT":
            key_list = key_words(index)
            while not file_clean[index + 1][0].startswith('*'):
                index += 1
                elements.append([int(num) for num in file_clean[index]])
        elif first == "*EL FILE":
            key_list = key_words(index)
            index = data_lines(index)
        elif first == "*ELSET":
            key_list = key_words(index)
            name = key_list[0]
            while not file_clean[index + 1][0].startswith('*'):
                index += 1
                if file_clean[index][0] == "EVOLUMES":
                    set = "EVOLUMES"
                else:
                    set = int(file_clean[index][0])
                elset.append([set, name.replace("ELSET=", "")])

        elif first == "*SURFACE":
            index = data_lines(index)
        elif first == "*DLOAD":
            while not file_clean[index + 1][0].startswith('*'):
                index += 1
                a = file_clean[index]
                sload.append([int(a[0]), a[1], int(a[2])])
        elif first == "*CONTACT":
            index = data_lines(index)
        elif first == "*FRICTION":
            index = data_lines(index)
        elif first == "*MPC":
            index = data_lines(index)
        elif first == "*SOLID":
            index = data_lines(index)
        elif first == "*SOLID SECTION":
            index = data_lines(index)
        elif first == "*ORIENTATION":
            index = data_lines(index)
        elif first == "*MATERIAL":
            key_list = key_words(index)
            name = key_list[0]
            index+=1
            if file_clean[index][0] == "*ELASTIC":
                index+=1
                a = file_clean[index]
                matset.append([name.replace("NAME=", ""),[float(a[0]), float(a[1])]])




        elif first == "*ELASTIC":
            index = data_lines(index)
        elif first == "*NSET":
            key_list = key_words(index)
            name = key_list[0]
            while not file_clean[index + 1][0].startswith('*'):
                index += 1
                nset.append([int(file_clean[index][0]), name.replace("NSET=", "")])
        elif first == "*BOUNDARY":
            key_list = key_words(index)
            while not file_clean[index + 1][0].startswith('*'):
                index += 1
                a = file_clean[index]
                bcon.append([a[0],int(a[1])])


        elif first == "*PRE-TENSION":
            index = data_lines(index)
        elif first == "*STEP":
            index = data_lines(index)
        elif first == "*STATIC":
            index = data_lines(index)
        elif first == "*CLOAD":
            index = data_lines(index)
        elif first == "*PRINT":
            index = data_lines(index)
        elif first == "*NODE PRINT":
            index = data_lines(index)
        elif first == "*END STEP":
            index = data_lines(index)
        elif first == "*END":
            pass
        elif first.startswith('**'):
            pass
        elif first == "":
            print("remove space between * and first")
        else:
            print("UNKNOWN")

print("_____________________________________________________")
print(nodes)
print(elements)
print(elset)
print(nset)
print(sload)
print(bcon)
print(matset)