#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a GUI that uses PyQt6 (version 6.3.0). 

It will allow the user to send a TTL pulse of a user-defined width
to any of the digital outputs (SMAs or parallel port). 

Created on Thu May  5 14:56:04 2022 

@author: kyleloizos
"""

import sys
import xipppy as xp
from time import sleep
from PyQt6.QtGui import QRegularExpressionValidator
from PyQt6.QtCore import QRegularExpression
from PyQt6.QtWidgets import ( 
    QApplication, 
    QWidget, 
    QPushButton,
    QCheckBox,
    QLineEdit,
    QLabel,
    QFormLayout,
    QVBoxLayout,
    QHBoxLayout
)

class MyApp(QWidget):
    def __init__(self):
        super().__init__()
        self.localParams()
        self.initUI()
        
    def localParams(self):
        self.connected = 0
        
    def initUI(self):
        # Set window properties
        self.setWindowTitle('Digital Output')
        self.resize(300, 250)
        
        # Create top level layout
        layout = QVBoxLayout()
        self.setLayout(layout)
        
        # Define button to connect/disconnect from processor
        self.connect_button = QPushButton("Connect to processor")
        self.connect_button.clicked.connect(self.connectToProcessor)
        
        # Define error label if xipppy can't connect to processor
        self.connectionError()
        
        # Define form for setting pulse width        
        pulsewidth_layout = QFormLayout()
        self.pulsewidth_text = QLineEdit("0.001")
        pulsewidth_layout.addRow("TTL width (s):", self.pulsewidth_text)
        
        # Define group of checkboxes for choosing digital output(s)
        digouts_layout = QVBoxLayout()
        self.digouts = [QCheckBox("SMA 1"), QCheckBox("SMA 2"), QCheckBox("SMA 3"), 
                        QCheckBox("SMA 4")]
        
        self.digouts_layout = QVBoxLayout()
        for index, elem in enumerate(self.digouts):
            digouts_layout.addWidget(self.digouts[index])
        
        # Define checkbox to enable parallel port input & text entry form
        parallel_layout = QHBoxLayout()
        self.parallel_enable = QCheckBox("Parallel port")
        parallel_layout.addWidget(self.parallel_enable)
        parallel_layout.addLayout(self.defineHexTextEntry())

        # Define button to send TTL
        self.sendTTL_button = QPushButton("Send TTL")
        self.sendTTL_button.clicked.connect(self.sendTTL)
        self.sendTTL_button.setEnabled(False)
        
        # Assemble window
        layout.addStretch()
        layout.addWidget(self.connect_button)
        layout.addWidget(self.connect_error)
        layout.addLayout(pulsewidth_layout)
        layout.addLayout(digouts_layout)
        layout.addLayout(parallel_layout)
        layout.addWidget(self.sendTTL_button)
        layout.addStretch()
        
    def connectionError(self):
        self.connect_error = QLabel("Connection failed. Check power and connection.\n")
        self.connect_error.setStyleSheet("color: red")
        self.connect_error.setHidden(True)
    
    def connectToProcessor(self):
        # Close connection
        self.connect_error.setHidden(True)
        xp._close()
        sleep(0.001)
        
        # Connect to processor (Try UDP then TCP). 
        # Enable TTL button when connected
        if self.connected==0:
            try:
                xp._open()
            except:
                try:
                    xp._open(use_tcp=True)
                except:
                    self.connected = 0
                    self.connect_error.setHidden(False)
                    print("Failed to connect to processor.")
                else:
                    self.connected = 1
                    print("Connected over TCP")
                    self.sendTTL_button.setEnabled(True)
                    self.connect_button.setText("Disconnect")
            else:
                self.connected = 1
                print("Connected over UDP")
                self.sendTTL_button.setEnabled(True)
                self.connect_button.setText("Disconnect")
        
        # If already connected, change button from 'connect' to 'disconnect'
        else:
            print("Closed connection")
            xp._close()
            self.sendTTL_button.setEnabled(False)
            self.connect_button.setText("Connect to processor")
            self.connected = 0
        sleep(0.001)
            
    def sendTTL(self):        
        digout_cmd = "xp.digout("
        out_channels = ""
        out_output = ""
        out_output2 = ""
        
        # Find all selected channels, populate digout subfields in strings
        # Start with parallel port
        if self.parallel_enable.isChecked():
            out_channels = "4"
            out_output = str(int(self.parallel_text.text(),16))
            out_output2 = "0"
        
        # Check all SMAs
        for index, elem in enumerate(self.digouts):
            if elem.isChecked():
                if out_channels != "":
                    out_channels = out_channels + "," + str(index)
                    out_output = out_output + "," + "1"
                    out_output2 = out_output2 + "," + "0"
                else:    
                    out_channels = out_channels + str(index)
                    out_output = out_output + "1"
                    out_output2 = out_output2 + "0"
             
        # Build digout command strings         
        digout_cmd1 = digout_cmd + "[" + out_channels + "], [" + out_output + "])"
        digout_cmd2 = digout_cmd + "[" + out_channels + "], [" + out_output2 + "])"
        
        # Execute TTL
        exec(digout_cmd1)
        sleep(float(self.pulsewidth_text.text()))
        exec(digout_cmd2)
        
    def defineHexTextEntry(self):
        parallel_form = QFormLayout()
        self.parallel_text = QLineEdit("CAFE")
        validator = QRegularExpressionValidator(QRegularExpression("[0-9A-Fa-f]{1,4}"))
        self.parallel_text.setValidator(validator)
        parallel_form.addRow(self.parallel_text) 
        
        # Disable form if checkbox is unchecked
        self.parallel_text.setEnabled(False)
        self.parallel_enable.toggled.connect(self.parallel_checked)
        
        return parallel_form    
    
    def parallel_checked(self):
        # If parallel port checkbox is checked, enable text entry
        if self.parallel_enable.isChecked():
            self.parallel_text.setEnabled(True)
        else:
            self.parallel_text.setEnabled(False)

if __name__ == "__main__":
    xp._close()
    app = QApplication([])

    window = MyApp()
    window.show()

    sys.exit(app.exec())