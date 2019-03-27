import sqlite3
import re


class DatabaseManager:

    def __init__(self, database):
        self.coms = {"su": "surfactant",
                "se": "sensitizer",
                "p": "photolysis",
                "d": "measured_time",
                "n": "name"}
        self.yyyy = re.compile("^\d{4}$")
        self.yyyymm = re.compile("^\d{6}$")
        self.yyyymmdd = re.compile("^\d{8}$")

        self.db = sqlite3.connect(database)
        self.cur = self.db.cursor()

        self.current_request = None


    def do_command(self, command):
        self.cur.execute(command)

    def translate_command(self, command, view="def", table="regular_sfg", base_table="sfg"):

        out = None
        temp = command.split(" ")
        com = temp[0]

        if com not in self.coms:
            raise ValueError(f'Invalid command "{com}"!')

        if com in ("su", "se", "name"):
            out = f'SELECT * FROM {table} WHERE {self.coms[com]}'
            if len(temp) == 2:
                out += f'="{temp[1]}"'

            else:
                out += ' in ('
                for element in temp[1:-1]:
                    out += f'"{element}",'

                out += f'"{temp[-1]}")'

        elif com == "p":

            if len(temp) > 2 or temp[1] not in ("y", "n"):
                raise ValueError("Invalid request for photolized samples!")
            else:
                out = f'SELECT * FROM {table} WHERE {self.coms[com]} '
                if temp[1] == "y":
                    out += 'NOT NULL'
                else:
                    out += 'IS NULL'

        elif com == "d":

            if 2 < len(temp) > 3:
                raise ValueError("Invalid date request!")

            out = f"""SELECT * FROM {table} INNER JOIN {base_table} ON {base_table}.id={table}.specid
            WHERE {base_table}.{self.coms[com]} BETWEEN """
            date1 = self.process_date(temp[1])

            if len(temp) == 3:
                date2 = self.process_date(temp[2])

                if date1 and date2:
                    out += date1 + " AND " + date2

                else:
                    raise ValueError(f'Invalid date format in "{temp[1]}" or \
                                                                "{temp[2]}"')

            elif len(temp) == 2:
                if date1:
                    out += date1 + " AND " + date1
                else:
                    raise ValueError(f'"Invalid date format "{temp[1]}"')

        if view is not None:
            out = f'CREATE TEMP VIEW {view} AS ' + out

        print(out)
        return out

    def construct_objects(self, records):
        pass

    def process_date(self, date):
        out = None

        if re.match(self.yyyy, date):
            out = f'"{date}-01-01"'

        elif re.match(self.yyyymm, date):
            out = f'"{date[0:4]}-{date[4:6]}-01"'

        elif re.match(self.yyyymmdd, date):
            out = f'"{date[0:4]}-{date[4:6]}-{date[6:]}"'

        return out

    def show_view(self):
        return self.current_request

    def new_request(self, command, view="def", table="regular_sfg"):
        command = self.translate_command(command, view=view, table=table)
        self.do_command(command)

        if view is not None:
            command2 = f' SELECT id,name FROM {view}'
            self.do_command(command2)
            self.current_request = self.cur.fetchall()

        else:
            self.current_request = self.cur.fetchall()

    def restore(self, view):
        command = f'SELECT id,name FROM {view}'
        self.do_command(command)
        self.current_request = self.cur.fetchall()



class SessionControlManager:

    def __init__(self):
        self.dm = DatabaseManager("test.db")
        self.active_view = None
        self.views = []
        self.refined = False
        self.anonymous_count = 0

    def get(self, command, view="def"):
        self.dm.new_request(command, view=view)
        self.active_view = view
        self.refined = False

    def ref(self, command, view=None):
        if view is None:
            self.anonymous_count += 1
            view = "unnamed" + str(self.anonymous_count)

        self.dm.new_request(command, table=self.active_view, view=view)
        self.refined = True

        if view is not None:
            self.views.append(self.active_view)
            self.active_view = view

    def rec(self, view="default"):
        if view == "default":
            self.dm.restore(self.views[-1])
            self.active_view = self.views [-1]
            self.views.remove(self.views[-1])
        else:
            for item in self.views:
                if item == view:
                    self.dm.restore(item)
                    self.active_view = item
                    self.views.remove(item)

    def rm(self, command):
        pass

    def keep(self, command):
        pass

    def show(self):
        temp = self.dm.show_view()
        print(f'Current view of name {self.active_view} with refinement state {self.refined}: \n')

        for item in temp:
            print(f'id: {item[0]}, name: {item[1]}')

    def plot(self):
        pass

    def clear(self):
        self.active_view = None
        self.refined = False
        self.dm.current_request = []

#testcode section

S = SessionControlManager()
S.get("su BX12")
S.ref("p y", view="photo")
S.ref("d 20171101 20171231")
S.show()
S.rec("def")
print(S.views, S.active_view)
