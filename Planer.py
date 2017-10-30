class Planer:

    """A class to perform task management and plan the research tasks for a period of measurements"""

    def __init__(self):

        self.tasks = []
        self.count = 0
        self.refresh()
        print("Planer initialized. The current task count is "+str(self.count))
        self.show_tasks()

    def show_tasks(self):

        index = 1
        for i in self.tasks:
            print(str(index)+" : "+i)
            index += 1

    def add_task(self, string):

        with open("tasks", "a") as infile:
            infile.write(string+"\n")
        self.refresh()

    def refresh(self):
        self.tasks = []
        with open("tasks", "r") as infile:

            for line in infile:
                self.tasks.append(line)
                self.count = len(self.tasks)

    def done(self, number):

        print(str(self.tasks[number-1])+"was removed!")
        del self.tasks[number-1]
        self.update_file()

    def update_file(self):

        with open("tasks", "w") as outfile:
            for task in self.tasks:
                outfile.write(task)



#test code section

P = Planer()
#P.add_task("task added")
P.show_tasks()
P.done(4)
P.show_tasks()