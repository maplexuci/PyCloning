from Bio.Seq import Seq
from tkinter import *
from tkinter import scrolledtext
from PIL import ImageTk, Image
# from tkinter.font import Font


class main:
    """The main app class"""
    def __init__(self):
        root = Tk()
        root_window = Root(root)
        return None


class Root:
    """Root window class including the widgets and the functions
    for responding the widgets.
    """
    def __init__(self, root):
        # Root window configration
        self.root = root
        self.root.title("PyCloning")
        self.root.geometry("700x400")
        self.root.resizable(width=False, height=False)
        self.root.configure(bg="#ebf5fc")
        self.root.columnconfigure(0, weight=2)
        self.root.columnconfigure((1, 2), weight=1)

        # Call functions to GUI interface
        self._create_menu()
        self._create_display()
        self._create_button_and_label()
        self._create_statusBar()

        # Start mainloop
        self.root.mainloop()

    def _create_menu(self):
        # Create menubar first
        self.menubar = Menu(self.root)
        self.root.config(menu=self.menubar)

        # Create 'File' menu inside the menubar
        # 'tearoff=False' ensures the menu is not a floatting menu
        self.file_menu = Menu(self.menubar, tearoff=False)
        # Adding cascade to 'File' menu.
        self.menubar.add_cascade(label="File", menu=self.file_menu)
        self.file_menu.add_command(label="New DNA File...", command=lambda: NewDNA(self))  # 'self' is passed to the 'parent' argument in class NewDNA constructor.
        self.file_menu.add_command(label="New Protein File...")
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Open Files...")
        self.file_menu.add_command(label="Open Recent File")
        self.file_menu.add_command(label="Close", command=quit)
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Import")
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Exit", command=quit)

        # Create 'Cloning' menu
        self.cloning_menu = Menu(self.menubar, tearoff=False)
        self.menubar.add_cascade(label="Cloning", menu=self.cloning_menu)
        self.cloning_menu.add_command(label="PCR...")
        self.cloning_menu.add_command(label="Overlapping PCR")
        self.cloning_menu.add_command(label="Mutagenesis...")
        self.cloning_menu.add_separator()
        self.cloning_menu.add_command(label="Restriction Enzyme Cloning")
        self.cloning_menu.add_command(label="Linear Ligation")
        self.cloning_menu.add_separator()
        self.cloning_menu.add_command(label="TA or GC Cloning")
        self.cloning_menu.add_command(label="TOPO® Cloning")
        self.cloning_menu.add_command(label="Gateway® Cloning")
        self.cloning_menu.add_command(label="Gibson Assembly®")
        self.cloning_menu.add_command(label="In-Fusion® Cloning")
        self.cloning_menu.add_command(label="NEBuilder® HiFi DNA Assembly")

        # Create 'Reference' menu
        self.ref_menu = Menu(self.menubar, tearoff=False)
        self.menubar.add_cascade(label="Reference", menu=self.ref_menu)
        self.ref_menu.add_command(label="Restriction Enzymes")
        self.ref_menu.add_command(label="Enzyme Database")
        self.ref_menu.add_separator()
        self.ref_menu.add_command(label="Common Features")
        self.ref_menu.add_command(label="Letter Codes")
        self.ref_menu.add_separator()
        self.ref_menu.add_command(label="Genetic Code Tables")
        self.ref_menu.add_command(label="Codon Usage Tables")

    def _create_display(self):
        self.VERSION = "V0.1"
        self.ver_title = Label(self.root, text="PyCloning " + self.VERSION, font=("Courier", 24, "bold"), bg="#ebf5fc")
        self.ver_title.grid(row=0, column=0, columnspan=3, padx=220, pady=10, sticky=E+W)

        self.canvas = Canvas(self.root, bg="#ebf5fc", highlightthickness=0, width=260, height=260)
        self.canvas.grid(row=1, column=0, rowspan=5, padx=3, pady=9)
        self.logo = ImageTk.PhotoImage(Image.open("imgs/logo.png"))
        self.canvas.create_image(130, 130, image=self.logo)

    def _create_button_and_label(self):
        self.ndf_img = PhotoImage(file="imgs/dna_icon_30.png")
        self.npf_img = PhotoImage(file="imgs/protein_icon_30.png")
        self.open_img = PhotoImage(file="imgs/folder_empty_icon_30.png")
        self.orf_img = PhotoImage(file="imgs/folder_file_icon_30.png")
        self.import_img = PhotoImage(file="imgs/import_icon_30.png")
        self.btn_ndf = Button(self.root, image=self.ndf_img, width=35, height=35)
        self.btn_npf = Button(self.root, image=self.npf_img, width=35, height=35)
        self.btn_open = Button(self.root, image=self.open_img, width=35, height=35)
        self.btn_orf = Button(self.root, image=self.orf_img, width=35, height=35)
        self.btn_import = Button(self.root, image=self.import_img, width=35, height=35)

        self.lbl_ndf = Label(self.root, text="New DNA File...", bg="#ebf5fc")
        self.lbl_npf = Label(self.root, text="New Protein File...", bg="#ebf5fc")
        self.lbl_open = Label(self.root, text="Open", bg="#ebf5fc")
        self.lbl_orf = Label(self.root, text="Open Recent File", bg="#ebf5fc")
        self.lbl_import = Label(self.root, text="Import", bg="#ebf5fc")

        self.btn_ndf.grid(row=1, column=1, pady=5)
        self.btn_npf.grid(row=2, column=1, pady=5)
        self.btn_open.grid(row=3, column=1, pady=5)
        self.btn_orf.grid(row=4, column=1, pady=5)
        self.btn_import.grid(row=5, column=1, pady=5)

        self.lbl_ndf.grid(row=1, column=2, pady=10, sticky=W)
        self.lbl_npf.grid(row=2, column=2, pady=10, sticky=W)
        self.lbl_open.grid(row=3, column=2, pady=10, sticky=W)
        self.lbl_orf.grid(row=4, column=2, pady=10, sticky=W)
        self.lbl_import.grid(row=5, column=2, pady=10, sticky=W)

    def _create_statusBar(self):
        self.status = Label(self.root, bd=1, relief=SUNKEN, height=3)
        self.status.grid(row=6, column=0, columnspan=3, sticky=W+E)

    def hide(self):
        """Hide the root window."""
        self.root.withdraw()

    def show(self):
        """Show the root window from the hide status"""
        self.root.update()
        self.root.deiconify()

    def onClosing(self, window):
        """Respond to the toplevle window closing event:
        Close the current window and show the root window.
        """
        window.destroy()
        self.show()


class NewDNA:
    def __init__(self, parent):
        """Create new toplevel window for new DNA file
        Variable:
            parent: the widget class, from which this NewDNA class is called.

        Returns:
            [type]: [None]
        """
        # 'parent' refers to the Root class instance.
        # 'self.parent=parent' associates the Root instance with NewDNA instance.
        self.parent = parent
        self.parent.hide()  # This will call the hide() method in Root class.
        self.dna_window = Toplevel()
        self.dna_window.title("New DNA File")
        self.dna_window.geometry("500x400")

        # 'grab_set()' and `focus()` ensure the toplevel is active (focused).
        self.dna_window.grab_set()
        self.dna_window.focus()

        # Call functions to build GUI interface
        self._create_infoLabel()
        self._create_seqInput()
        self._create_seqLen_display()
        self._create_fileName()
        self._create_buttons()

    def _create_infoLabel(self):
        lbl = Label(self.dna_window, text="Input your sequence:")
        lbl.pack(padx=10, pady=(10, 0), anchor=W)

    def _create_seqInput(self):
        # Create a scrolledtext widget.
        self.new_dna_seq = scrolledtext.ScrolledText(
                                self.dna_window, wrap=WORD,
                                width=80, height=15, font=("Courier", 11),
                                )

        # 'expand=True' and 'fill=BOTH' ensure that
        # the Text widget change size along with window resizing.
        self.new_dna_seq.pack(padx=10, expand=True, fill=BOTH, anchor=W)
        self.new_dna_seq.focus_set()

        # Bind the scrooedtext widget with 'self.changed' function
        # to keep track if the content in the widget is 'Modified'.
        self.new_dna_seq.bind('<<Modified>>', self.changed)
        self.new_dna_seq.bind('<KeyRelease>', self.onValidate)

    def _create_seqLen_display(self):
        # Add another Label widget to display the sequence length.
        self.seq_len_lbl = Label(self.dna_window)
        self.seq_len_lbl.pack(padx=10, anchor=W)

    def _create_fileName(self):
        # Add Label and Entry widget in a Frame for entering filename.
        # First, create Frame widget
        filename_Frame = Frame(self.dna_window)

        # Create text Label
        Label(filename_Frame, text="File Name: ").pack(side=LEFT)

        # Define and use variable type for Entry widget.
        self.name_var = StringVar()
        seq_filename = Entry(filename_Frame, textvariable=self.name_var)
        seq_filename.pack(side=LEFT)

        # Pack the Frame in the end after all widgets in Frame are packed.
        filename_Frame.pack(padx=10, anchor=W)

    def _create_buttons(self):
        # Respond to the 'Cancel' button.
        btn_cancel = Button(self.dna_window, text="Cancel", width=10,
                            command=lambda: self.parent.onClosing(self.dna_window))
        btn_cancel.pack(padx=10, pady=10, side=RIGHT)

        # Add 'OK' button to read sequence
        self.btn_ok = Button(self.dna_window, text="OK", width=10, state=DISABLED,
                             command=self.readSeq)
        self.btn_ok.pack(padx=10, pady=10, side=RIGHT)

        #  Respond to 'close window' event
        self.dna_window.protocol("WM_DELETE_WINDOW",
                                 lambda: self.parent.onClosing(self.dna_window))

        return None

    def changed(self, event):
        """Function to keep track the changes in Text and reflect the changes in the Label"""
        # Only display the length of sequence when it is not empty.
        text = ''

        # '1.0' means the first row (number '1'), the first column (index '0').
        text_len = len(self.new_dna_seq.get(1.0, END).rstrip())

        # Give a 'flag' to track if the ScrolledText widget
        # is modified using '.edit_modified()' function.
        # flag has two values, 1 or Trun means modified;
        # 0 or False means unmodified.
        flag = self.new_dna_seq.edit_modified()
        if flag == 1:
            if text_len != 0:
                text = f"{str(text_len)} bp"
            else:
                text = ''

            # Display the sequence length on the Label.
            # We need to display it first before reset
            # the modification flag to 0.
            # otherwise, it won't display, because the reset statement
            # immediately call the changed() function and the condition in
            # if statement will not meet.
            self.seq_len_lbl.config(text=text)

            # Reset the modified flag to False (0), for detecting
            # new modification. Note that, this also means the modification
            # state is changes again, so, this will call the changed()
            # function once again. How ever, we set a control condition
            # 'if self.flag == 1', this will ensure the code inside of this
            # contition statement not excecute again.
            self.new_dna_seq.edit_modified(False)        

    def onValidate(self, event):
        char = ('A', 'a', 'T', 't', 'G', 'g', 'C', 'c', 'W', 'w', 'S', 's',
                'M', 'm', 'K', 'k', 'R', 'r', 'Y', 'y', 'B', 'b', 'D', 'd',
                'H', 'h' 'V', 'v', 'N', 'n')
        seq = self.new_dna_seq.get(1.0, END).rstrip()
        if len(seq) != 0:
            self.btn_ok.config(state=NORMAL)
        else:
            self.btn_ok.config(state=DISABLED)
        for letter in seq:
            if letter not in char:
                self.new_dna_seq.delete(1.0, END)
                self.new_dna_seq.insert(1.0, seq.replace(letter, ''))

                # Update seq to ensure loops through the whole sequence.
                seq = self.new_dna_seq.get(1.0, END).rstrip()

    def get_Filename(self):
        """Get the content in an Entry as string.

        Args:
            name_var (var): The variable type defined for the Entry.

        Returns:
            Str: The content in an Entry widget.
        """

        return self.name_var.get()

    def readSeq(self):
        self.dnaSeq = Seq(self.new_dna_seq.get(1.0, "end-1c"))
        self.dna_window.destroy()
        workWindow = WorkingWindow(self)

    def call_root_function(self, window):
        """Call onClosing() function in with a Root instance"""
        self.parent.onClosing(window)


class WorkingWindow:
    """Class for the main working window, where you manipulate sequences.
    """
    def __init__(self, parent):
        """Create the working window GUI interface.
        """
        self.parent = parent
        self.seq_window = Toplevel()
        self.title = parent.get_Filename()
        self.seq_window.title(self.title)
        self.seq_window.state('zoomed')  # Make the window maximized.

        # To call onClosing function() in Root, we first call
        # 'call_root_function()' in NewDNA class using 'self.parent', which refers to NewDNA class instance.
        # Then in 'call_root_function()', we call 'onClosing()' in Root, using 'self.parent', which refers to Root class instance.
        self.seq_window.protocol("WM_DELETE_WINDOW",
                                 lambda: self.parent.call_root_function(self.seq_window))

        # Call functions to build GUI interface
        self._seqEditor()
        self._workWindowMenu()

    def _workWindowMenu(self):
        # Create the main menubar first
        self.menubar = Menu(self.seq_window)
        self.seq_window.config(menu=self.menubar)

        # Create items for 'File' menu.
        # 'tearoff=False' ensures the menu is not a floatting menu
        self.file_menu = Menu(self.menubar, tearoff=False)
        self.menubar.add_cascade(label="File", menu=self.file_menu)
        self.file_menu.add_command(label="New DNA File...")
        self.file_menu.add_command(label="New Protein File...")
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Open Files...")
        self.file_menu.add_command(label="Open Recent File")
        self.file_menu.add_command(label="Close", command=quit)
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Import")
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Exit", command=quit)

        # Create items for 'Edit' menu.
        self.edit_menu = Menu(self.menubar, tearoff=False)
        self.menubar.add_cascade(label="Edit", menu=self.edit_menu)
        self.edit_menu.add_command(label="Undo")
        self.edit_menu.add_command(label="Redo")
        self.edit_menu.add_separator()
        self.edit_menu.add_command(label="Cut")
        self.edit_menu.add_command(label="Copy")
        self.edit_menu.add_command(label="Paste")
        self.edit_menu.add_command(label="Delete")
        self.edit_menu.add_separator()
        self.edit_menu.add_command(label="Select All")
        self.edit_menu.add_command(label="Select Range...")
        self.edit_menu.add_separator()
        self.edit_menu.add_command(label="Make Uppercase")
        self.edit_menu.add_command(label="Make Lowercase")
        self.edit_menu.add_command(label="Set Color...")
        self.edit_menu.add_separator()
        self.insert_sub = Menu(self.edit_menu, tearoff=False)
        self.edit_menu.add_cascade(label="Insert", menu=self.insert_sub)
        self.insert_sub.add_command(label="Base...")
        self.insert_sub.add_command(label="Codon...")
        self.insert_sub.add_command(label="Restriction Site...")
        self.insert_sub.add_command(label="Feature...")
        self.edit_menu.add_separator()
        self.find_sub = Menu(self.edit_menu, tearoff=False)
        self.edit_menu.add_cascade(label="Find", menu=self.find_sub)
        self.find_sub.add_command(label="Find DNA Sequence")
        self.find_sub.add_command(label="Find Protein Sequence")
        self.find_sub.add_command(label="Find Enzyme / Fearure / Primer")
        self.edit_menu.add_command(label="Go To...")
        self.edit_menu.add_separator()
        self.edit_menu.add_command(label="Preference")
        self.language_sub = Menu(self.edit_menu, tearoff=False)
        self.edit_menu.add_cascade(label="Language", menu=self.language_sub)
        self.language_sub.add_command(label="English")
        self.language_sub.add_command(label="Chinese")

        # Create items for 'View menu
        self.view_menu = Menu(self.menubar, tearoff=False)
        self.menubar.add_cascade(label="View", menu=self.view_menu)

        # Create items for 'Enzymes' menu
        self.enzyme_menu = Menu(self.menubar, tearoff=False)
        self.menubar.add_cascade(label="Enzymes", menu=self.enzyme_menu)
        self.enzyme_menu.add_command(label="Choose Enzymes...")
        self.enzyme_menu.add_command(label="Show All Enzyme Sites")
        self.enzyme_menu.add_command(label="Show Unique Enzyme Sites")
        self.enzyme_menu.add_command(label="Hide All Enzymes")
        self.enzyme_menu.add_command(label="Noncutters...")
        self.enzyme_menu.add_separator()
        self.enzyme_menu.add_command(label="Restriction Enzymes...")
        self.enzyme_menu.add_command(label="Enzyme Database...")

        # Create items for 'Features' menu
        self.feature_menu = Menu(self.menubar, tearoff=False)
        self.menubar.add_cascade(label="Features", menu=self.feature_menu)
        self.feature_menu.add_command(label="Add Feature...")
        self.feature_menu.add_command(label="Edit Feature...")
        self.feature_menu.add_command(label="Duplicate Feature...")
        self.feature_menu.add_command(label="Remove Feature")
        self.feature_menu.add_separator()
        self.feature_menu.add_command(label="Create Feature Segment...")
        self.feature_menu.add_command(label="Delete Feature Segment...")
        self.feature_menu.add_command(label="Merge Feature Segments...")
        self.feature_menu.add_separator()
        self.feature_menu.add_command(label="Feature Color...")
        self.feature_menu.add_command(label="Feature List...")
        self.feature_showhide_sub = Menu(self.feature_menu, tearoff=False)
        self.feature_menu.add_cascade(label="Show/Hide Features", menu=self.feature_showhide_sub)
        self.feature_showhide_sub.add_command(label="Show All Features")
        self.feature_showhide_sub.add_command(label="Show Selected Features")
        self.feature_showhide_sub.add_command(label="Hide All Features")
        self.feature_showhide_sub.add_command(label="Hide Selected Features")
        self.feature_menu.add_separator()
        self.importFeature_sub = Menu(self.feature_menu, tearoff=False)
        self.feature_menu.add_cascade(label="Import Features", menu=self.importFeature_sub)
        self.importFeature_sub.add_command(label="Import Features from a SnapGene File...")
        self.importFeature_sub.add_command(label="Import Features from a BED File...")
        self.importFeature_sub.add_command(label="Import Features from a GFF3 File...")
        self.importFeature_sub.add_command(label="Import Features from a GTF File...")
        self.feature_menu.add_command(label="Export Feature Data...")

        # Create items in 'Primer' menu
        self.primer_menu = Menu(self.menubar, tearoff=False)
        self.menubar.add_cascade(label="Primers", menu=self.primer_menu)
        self.primer_menu.add_command(label="Add Primer...")
        self.primer_menu.add_command(label="Edit Primer...")
        self.primer_menu.add_command(label="Duplicate Primer...")
        self.primer_menu.add_command(label="Remove Primer")
        self.primer_menu.add_separator()
        self.primer_menu.add_command(label="Primer Color...")
        self.primer_menu.add_command(label="Hybridization Parameters...")
        self.primer_menu.add_separator()
        self.primer_showhide_sub = Menu(self.primer_menu, tearoff=False)
        self.primer_menu.add_cascade(label="Show/Hide Primers", menu=self.primer_showhide_sub)
        self.primer_showhide_sub.add_command(label="Show All Primers")
        self.primer_showhide_sub.add_command(label="Show Selected primers")
        self.primer_showhide_sub.add_command(label="Hide All Primers")
        self.primer_showhide_sub.add_command(label="Hide Selected Primers")
        self.primer_menu.add_command(label="Primer List...")

        # Create items for 'Cloning' menu
        self.cloning_menu = Menu(self.menubar, tearoff=False)
        self.menubar.add_cascade(label="Cloning", menu=self.cloning_menu)
        self.cloning_menu.add_command(label="PCR...")
        self.cloning_menu.add_command(label="Overlapping PCR")
        self.cloning_menu.add_command(label="Mutagenesis...")
        self.cloning_menu.add_separator()
        self.cloning_menu.add_command(label="Restriction Enzyme Cloning")
        self.cloning_menu.add_command(label="Linear Ligation")
        self.cloning_menu.add_separator()
        self.cloning_menu.add_command(label="TA or GC Cloning")
        self.cloning_menu.add_command(label="TOPO® Cloning")
        self.cloning_menu.add_command(label="Gateway® Cloning")
        self.cloning_menu.add_command(label="Gibson Assembly®")
        self.cloning_menu.add_command(label="In-Fusion® Cloning")
        self.cloning_menu.add_command(label="NEBuilder® HiFi DNA Assembly")

        # Create items for 'Tools' menu
        self.tool_menu = Menu(self.menubar, tearoff=False)
        self.menubar.add_cascade(label="Tools", menu=self.tool_menu)

        # Create items for 'Window' menu
        self.window_menu = Menu(self.menubar, tearoff=False)
        self.menubar.add_cascade(label="Windows", menu=self.window_menu)

        # Create items for 'Help' menu
        self.help_menu = Menu(self.menubar, tearoff=False)
        self.menubar.add_cascade(label="Help", menu=self.help_menu)

    def _seqEditor(self):
        # Create a scrolledtext widget.
        self.dnaSeq = self.parent.dnaSeq

        self.seq = scrolledtext.ScrolledText(self.seq_window, wrap=WORD,
                                             font=("Consolas", 12), bg="#f5feff")
        # 'expand=True' and 'fill=BOTH' ensure that
        # the Text widget change size along with window resizing."
        self.seq.pack(padx=10, expand=True, fill=BOTH, anchor=W)

        row_num = len(self.dnaSeq)//150 + 1
        remain_seq_len = 0
        current_seq_len = 0
        for row in range(row_num):
            if row == 0:
                row_seq = self.dnaSeq[0 : (row+1)*150]
                current_seq_len = (row+1)*150 + remain_seq_len
            elif 0 < row <= row_num-2:
                row_seq = self.dnaSeq[row*150+1:(row+1)*150+1]
                current_seq_len = (row+1)*150 + remain_seq_len
            else:
                remain_seq_len = len(self.dnaSeq) % 150
                row_seq = self.dnaSeq[row*150+1:-1]
                current_seq_len = current_seq_len + remain_seq_len
            
            whitespace = ' '
            space = 150-len(row_seq)+8
            row_seq_com = self._complement(row_seq)
            row_display = f"{row_seq}{whitespace*space}{current_seq_len}\n{row_seq_com}\n"
            self.seq.insert(END, row_display)
            self.seq.insert(END, '\n')

            # self.seq.window_create(2.5, window=Button(self.seq, text="Feature"))  ## This is a way to add Sequence features. Need to seperate sequences
                                                                              ## into multi-lines for adding more widgets at desired index.

    def _complement(self, dna):
        # Return the complement sequence of the Seq object.
        self.dna = dna
        return self.dna.complement()
    
    def _reverseComplement(self, seq):
        # Return the reverse complement sequence of the Seq object.
        self.seq = seq
        return self.seq.reverse_complement()

    


if __name__ == '__main__':
    main()
