
class Person(object):
    """A person object"""

    def __init__(self, age, height, weight):
        self.age = age
        self.height = height
        self.weight = weight

    def how_old(self):
        print("He is {age} years old and he is {height} inches tall.".format(age=self.age, height=self.height))

    def bmi(self):
        bmivalue = float(self.weight)/float(self.height)**2 * 703
        print(bmivalue)

    def birth_year(self, year):
        self.birthyear = year - self.age
        print(self.birthyear)

class Teacher(Person):
    pass


paul = Teacher(47, 73, 160)

paul.how_old()

paul.bmi()

paul.birth_year(2014)