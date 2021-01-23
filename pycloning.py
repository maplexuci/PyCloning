from Bio.Seq import Seq
from tkinter import *
from PIL import ImageTk, Image


root = Tk()
root.title("PyCloning")
root.geometry("700x400")
root.resizable(width=False, height=False)
root.configure(bg="#ebf5fc")
root.columnconfigure(0, weight=2)
root.columnconfigure((1,2), weight=1)

# Create menu bar
menubar = Menu(root)
root.config(menu=menubar)

# Create 'File' menu
file_menu = Menu(menubar, tearoff=False)
menubar.add_cascade(label="File", menu=file_menu)
file_menu.add_command(label="New DNA File...")
file_menu.add_command(label="New Protein File...")
file_menu.add_separator()
file_menu.add_command(label="Open Files...")
file_menu.add_command(label="Open Recent File")
file_menu.add_command(label="Close", command=root.quit)
file_menu.add_separator()
file_menu.add_command(label="Import")
file_menu.add_separator()
file_menu.add_command(label="Exit", command=root.quit)

# Create 'Cloning' menu
cloning_menu = Menu(menubar, tearoff=False)
menubar.add_cascade(label="Cloning", menu=cloning_menu)
cloning_menu.add_command(label="PCR...")
cloning_menu.add_command(label="Overlapping PCR")
cloning_menu.add_command(label="Mutagenesis...")
cloning_menu.add_separator()
cloning_menu.add_command(label="Restriction Enzyme Cloning")
cloning_menu.add_command(label="Linear Ligation")
cloning_menu.add_separator()
cloning_menu.add_command(label="TA or GC Cloning")
cloning_menu.add_command(label="TOPO® Cloning")
cloning_menu.add_command(label="Gateway® Cloning")
cloning_menu.add_command(label="Gibson Assembly®")
cloning_menu.add_command(label="In-Fusion® Cloning")
cloning_menu.add_command(label="NEBuilder® HiFi DNA Assembly")

# Create 'Reference' menu
ref_menu = Menu(menubar, tearoff=False)
menubar.add_cascade(label="Reference", menu=ref_menu)
ref_menu.add_command(label="Restriction Enzymes")
ref_menu.add_command(label="Enzyme Database")
ref_menu.add_separator()
ref_menu.add_command(label="Common Features")
ref_menu.add_command(label="Letter Codes")
ref_menu.add_separator()
ref_menu.add_command(label="Genetic Code Tables")
ref_menu.add_command(label="Codon Usage Tables")


VERSION = "V0.1"
ver_title = Label(root, text="PyCloning " + VERSION, font=("Courier", 24, "bold"), bg="#ebf5fc")
ver_title.grid(row=0, column=0, columnspan=3, padx=220, pady=10, sticky=E+W)

logo = ImageTk.PhotoImage(Image.open("imgs/logo.jpg"))
logo_label = Label(image=logo)
logo_label.grid(row=1, column=0, rowspan=5, padx=3, pady=9)

status = Label(root, bd=1, relief=SUNKEN, height=3)
status.grid(row=6, column=0, columnspan=3, sticky=W+E)

pixelVirtual = PhotoImage(width=1, height=1)
btn_ndf = Button(root, image=pixelVirtual, width=35, height=35)
btn_npf = Button(root, image=pixelVirtual, width=35, height=35)
btn_open = Button(root, image=pixelVirtual, width=35, height=35)
btn_orf = Button(root, image=pixelVirtual, width=35, height=35)
btn_import = Button(root, image=pixelVirtual, width=35, height=35)

lbl_ndf = Label(root, text="New DNA File...", bg="#ebf5fc")
lbl_npf = Label(root, text="New Protein File...", bg="#ebf5fc")
lbl_open = Label(root, text="Open", bg="#ebf5fc")
lbl_orf = Label(root, text="Open Recent File", bg="#ebf5fc")
lbl_import = Label(root, text="Import", bg="#ebf5fc")

btn_ndf.grid(row=1, column=1, pady=5)
btn_npf.grid(row=2, column=1, pady=5)
btn_open.grid(row=3, column=1, pady=5)
btn_orf.grid(row=4, column=1, pady=5)
btn_import.grid(row=5, column=1, pady=5)

lbl_ndf.grid(row=1, column=2, pady=10, sticky=W)
lbl_npf.grid(row=2, column=2, pady=10, sticky=W)
lbl_open.grid(row=3, column=2, pady=10, sticky=W)
lbl_orf.grid(row=4, column=2, pady=10, sticky=W)
lbl_import.grid(row=5, column=2, pady=10, sticky=W)


root.mainloop()
