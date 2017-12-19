import time
import os


class LogManager:

    def __init__(self):
        self.day = time.strftime("%d")
        self.month = time.strftime("%m")
        self.year = time.strftime("%Y")
        self.timetag = self.year+self.month+self.day
        self.sep = "*"*45
        self.sep2 = "-"*45
        self.spectra = []

        #self.manageDirectory()
        print(self.timetag)
        with open(self.timetag+".log" ,"w") as outfile:

            outfile.write("Logfile from day:    "+self.timetag+"\n")
            outfile.write(self.sep+"\n")

    def newSample(self):
        for file in os.listdir(os.getcwd()):
            if file.endswith("py"):
                if file not in self.spectra:
                    self.spectra.append(file)
                    t = self.getTimestamp()+"\t\t"+file
                    self.writeToLog(t+"\n")
                    meta = self.askMeta()
                    self.writeToLog(meta[0]+"\n")
                    self.writeToLog(meta[1]+"\n")
                    self.writeToLog(self.sep2+"\n")

    def newComment(self, comment):
        self.writeToLog(self.sep+"\n")
        self.writeToLog(comment+"\n")
        self.writeToLog(self.sep+"\n")

    def countSpecs(selfs):
        pass

    def getTimestamp(self):
        t = time.strftime("%H"+":"+"%M")
        return t

    def manageDirectory(self):
        os.mkdir(self.timetag)
        os.chdir(self.timetag)

    def writeToLog(self, string):

        with open(self.timetag+".log", "a") as outfile:
            outfile.write(string)

    def askMeta(self):
        purpose = input("Enter, if necessary, a description of the measurement´s purpose: ")
        observations = input("Enter any important observations made during measurements: ")
        return (purpose, observations)

L = LogManager()
L.newSample()
