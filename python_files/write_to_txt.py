L = [1.23, 1.3432, 45.434, 5.645, 456.543]


f = open("f1.txt", "w")
print(L)

for ele in L:
    f.write(ele + "\n")

