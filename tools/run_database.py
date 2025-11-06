import csv

def get_run_values(run):
    csvfile =  open("tools/run_database.csv",newline='')
    reader = csv.DictReader(csvfile)
    str_run = str(run)
    found_run = False
    l_return = [-1,-1,-1,-1,-1]
    for line in reader:
        if line["run"] == str_run:
            l_return[0] = int(line["entries"])
            for i in range(2,6):
                val = line[reader.fieldnames[i]] # csv header
                l_return[i-1] = int(val) if val != '' else -1
            found_run = True
            break
    if found_run:
        return l_return
    else:
        print("ERROR: couldn't find the run in the database")
        return l_return


def get_run_entries(run):
    e, f, t, s, m = get_run_values(run)
    return e

def write_eventfilter(run,eventfilter):
    # load entire file into memory
    with open("tools/run_database.csv",newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        lines = list(reader)
        myheader = reader.fieldnames

    # find the line to change
    for line in lines:
        if line["run"] == str(run):
            line["eventfilter"] = eventfilter
            break

    # write the file from memory
    with open("tools/run_database.csv",'w',newline='') as csvfile:
        writer = csv.DictWriter(csvfile,fieldnames=myheader)
        writer.writeheader()
        writer.writerows(lines)
