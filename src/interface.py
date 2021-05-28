"""functions and modules to generate the GUI and interface with the rest of the code
"""

import tkinter as tk

if __name__ == "__main__":
    window = tk.Tk()
    greeting = tk.Label(text="Hello, Tkinter")
    greeting.pack()
    window.mainloop()
