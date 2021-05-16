# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 16:38:31 2021

@author: Main
"""

import random

class thing():
    def __init__(self, parent=None):
        self.value = random.randrange(1,100)
        self.parent = parent
        
    def remove_parent(self):
        """Connects thing to thinglist (parent)"""
        if self.parent != None:
            self.parent.remove_thing(self)
            self.parent = None
        
    def connect_parent(self, new_parent):
        """Connects thing to thinglist (parent)"""
        # Remove old parent first (if applicable).
        if self.parent != None:
            self.remove_parent()
        
        self.parent = new_parent
        self.parent.add_thing(self)

class thinglist():
    def __init__(self):
        self.thinglist = []
    
    def add_thing(self, thing):
        self.thinglist.append(thing)
        
    def remove_thing(self, thing):
        self.thinglist.remove(thing)

t = thinglist()
t1 = thing(parent=t)
t2 = thing(parent=t)


t.add_thing(t1)
t.add_thing(t2)

print(t.thinglist)
print(t1.parent)

# t1.remove_parent()

# print(t.thinglist)
# print(t1.parent)
