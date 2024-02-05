# list to store file lines
lines = []
# read file
for i in range(0,20): ##from 1 to number of files .log as index starts from 0
    with open(str(i+1)+".log", 'r') as fp: #dir/i+namefile.log
            # read an store all lines into list
            lines = fp.readlines()           
            #print("file", i+1, "read")
# Write file
    n = 152
    with open(str(i+1)+".log", 'w') as fp:#dir/i+namefile.data
        for number, line in enumerate(lines):
        # iterate each line
            if number >= n and number<n+1000: #from (first line of data)-1 to (last line of data)-1 as index starts from 0
                fp.write(line)
            if number == n+1000:
                fp.write(line.strip())
                        #print("Files ", i+1, "written")
    # with open(str(i+1)+".data", 'r') as fp:#dir/i+namefile.data
    #     lines = fp.readlines()
    #     del lines[-1]
    # with open(str(i+1)+".data", 'w') as fp:#dir/i+namefile.data
    #     for line in lines:
    #         fp.write(line)