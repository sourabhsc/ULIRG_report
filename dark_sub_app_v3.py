from tkinter import *
import tkinter
from tkinter import messagebox

import numpy as np
def plot_things(directory):
    return 0
def dark_subtraction(directory, input_hot, output, cent_x, cent_y, mask_radii ):
    cent_x =1
    cent_y =2

    sth = plot_things(directory)
    
    return np.random.random(1), sth
def process():

    try:
        directory =Entry.get(E1)
        input_hot =Entry.get(E2)
        output = Entry.get(E3)
        
        cent_x = Entry.get(E4)
        cent_y = Entry.get(E5)
        
        mask_radii = Entry.get(E6)
        print (type(directory))
        if  directory :
            #E10 = Entry(top, bd =5, bg = "green", fg = "blue")
            #E10.grid(row=10,column=1)
            #Entry.insert(E10,0,"enter the value of direcotyr")

            messagebox.showwarning("Warning","Please enter correct operator")


        print (input_hot, directory, output,cent_x , mask_radii )
        happy_factor, sth = dark_subtraction(directory, input_hot, output, cent_x, cent_y, mask_radii  )
        '''
        if happy_factor >0.5:
            E10 = Entry(top, bd =5, bg = "green", fg = "blue")
            E10.grid(row=10,column=1)
            Entry.insert(E10,0,"Done")
        else:
            
            E10 = Entry(top, bd =5, bg = "red", fg = "blue")
            E10.grid(row=10,column=1)
            Entry.insert(E10,0,"Stuck")
        '''
    except UnboundLocalError :

        messagebox.showwarning("Warning","Please enter correct operator")
    except ValueError :

        messagebox.showwarning("Warning","Please enter the value in positive integer")
def renew():
    E4.delete(0,END)
    
def clear():
    
    E1.delete(0,END)
    E2.delete(0,END)
    E3.delete(0,END)
    E4.delete(0,END)
    E5.delete(0,END)
    E6.delete(0,END)
    E10.delete(0,END)

def flash(event, main):
    main.after(100, lambda: Button.config(bg = 'lightgrey'))
def window(main):
    B=Button(main, bd = 10, bg = "cyan", text ="Submit",command= process).grid(row=12,column=1, sticky='se',)
    sth = process()
    Renew_button1 = Button(main, bg="cyan", bd = 10, text = 'New', command = renew).grid(row=12,column=2, sticky='se')
    Delete_button = Button(main,bg="red", bd = 10, text = 'Clear', command = clear).grid(row=12,column=4, sticky='se')
    main.bind( flash)

top = tkinter.Tk()
top.title("SBC Dark")
L1 = Label(top, text="SBC subtraction software",).grid(row=0,column=1, sticky='se',)
L2 = Label(top, text="Work dir",).grid(row=1,column=0, sticky='se',)
L3 = Label(top, text="Galaxy name",).grid(row=2,column=0, sticky='se',)
L4 = Label(top, text="Output file",).grid(row=3,column=0 ,sticky='se',)
L5 = Label(top, text="cent_x",).grid(row=4,column=0, sticky='se',)
L6 = Label(top, text="cent_y",).grid(row=5,column=0, sticky='se',)
L7 = Label(top, text="mask_radii",).grid(row=6,column=0, sticky='se',)
L10 = Label(top, text="Plotting?????",).grid(row=10,column=0, sticky='se',)

E1 = Entry(top, bd =5)
E1.grid(row=1,column=1)
E2 = Entry(top, bd =5)
E2.grid(row=2,column=1)
E3 = Entry(top, bd =5)
E3.grid(row=3,column=1)
E4 = Entry(top, bd =5)
E4.grid(row=4,column=1)
E5 = Entry(top, bd =5)
E5.grid(row=5,column=1)
E6 = Entry(top, bd =5)
E6.grid(row=6,column=1)

window(top)
top.mainloop()