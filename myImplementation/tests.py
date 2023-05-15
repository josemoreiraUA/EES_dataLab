arr = [1,2,3,4,5]

for i in range(len(arr)):
    print("-----------")
    for j in range(len(arr)):
        if j % 2 == 0:
            continue
        print("j",j)