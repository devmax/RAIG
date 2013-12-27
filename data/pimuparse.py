import csv


def parse(files):

    #header = ["time", "pimu_time", "roll", "pitch", "yaw1", "yaw2",
    #          "ax", "ay", "az"]
    obs = [list() for i in xrange(len(files))]

    for j in xrange(len(files)):
        print "Parsing file ", j+1
        with open(files[j]) as datafile:
            data = csv.reader(datafile)
            data.next()
            count = 0
            for row in data:
                obs[j].append([float(item) for item in row])
                count += 1
                if count % 10000 == 0:
                    print "\rRow: ", count,
    return obs

if __name__ == "__main__":
    files = ['pimu_2013-12-19_7_filtered.csv', 'pimu_2013-12-19_10_filtered.csv',
             'pimu_2013-12-19_13_filtered.csv']

    obs = parse(files)
