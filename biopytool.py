import PySimpleGUI as sg


sg.theme("LightGrey2")

# Create the window
main_window = sg.Window("BioOne", layout=[[]], size=(800, 600))

# Create an event loop
while True:
    event, values = main_window.read()
    # End program if user closes window or
    # presses the OK button
    if event == "OK" or event == sg.WIN_CLOSED:
        break

main_window.close()
