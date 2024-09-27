import PySimpleGUI as sg

sg.theme("DarkAmber") # Add a touch of color
#All the stuff inside your window
layout = [ [sg.Text("Some test on Row 1")], 
           [sg.Text("Enter some stuff here: "), sg.InputText()],
           [sg.Button("Done"), sg.Button("cancel")]]

#Create the window
window = sg.Window("Window Title", layout)
#Event Loop to process "events" and get the "values" of the inputs
while True:
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == 'cancel':
        break
    print("You entered ", values[0])

window.close()