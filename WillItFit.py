#!/usr/bin/env python

##       Will It Fit        ##
##                          ##
##     Copyright 2013,      ##
## University of Copenhagen ##
##                          ##
##  Martin Cramer Pedersen  ##
##       mcpe@nbi.dk        ##

## This file is part of WillItFit.                                      ##
##                                                                      ##
## WillItFit is free software: you can redistribute it and/or modify    ##
## it under the terms of the GNU General Public License as published by ##
## the Free Software Foundation, either version 3 of the License, or    ##
## (at your option) any later version.                                  ##
##                                                                      ##
## WillItFit is distributed in the hope that it will be useful,         ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of       ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        ##
## GNU General Public License for more details.                         ##
##                                                                      ##
## You should have received a copy of the GNU General Public License    ##
## along with WillItFit. If not, see <http://www.gnu.org/licenses/>.    ##

##  If you use this software in your work, please cite: ##
##                                                      ##
##  Pedersen, M. C., Arleth, L. & Mortensen, K.         ##
##  J. Appl. Cryst. 46(6), 1894-1898                    ##

## Libraries
import os
import re
import shutil
import sys
import subprocess

try:
    import wx
except:
    print ("")
    print ("***********************************************************************")
    print ("* WillItFit failed to import wxPython - is it correctly installed...? *")
    print ("***********************************************************************")
    print ("")
    print ("Press Enter...")
    raw_input()
    sys.exit(1)

import wx.lib.scrolledpanel

## Plotting
try:
    import matplotlib 
except:
    print ("")
    print ("*************************************************************************")
    print ("* WillItFit failed to import MatPlotLib - is it correctly installed...? *")
    print ("*************************************************************************")
    print ("")
    print ("Press Enter...")
    raw_input()
    sys.exit(1)

matplotlib.interactive(True)
matplotlib.use('WXAgg')

try:
    import pylab
except:
    print ("")
    print ("********************************************************************")
    print ("* WillItFit failed to import PyLab - is it correctly installed...? *")
    print ("********************************************************************")
    print ("")
    print ("Press Enter...")
    raw_input()
    sys.exit(1)

## Define main class and text
class MainCls(wx.Frame):
    def __init__(self, parent, id):
        # Store initial directory upon initialization
        self.InitialDirectoryStr = os.path.abspath(os.path.dirname(sys.argv[0]))
        os.chdir(self.InitialDirectoryStr)        
        
        # Initial widgets
        wx.Frame.__init__(self, parent, id, 'Will It Fit...?', size = (1000, 700))
     
        BoxSizer = wx.BoxSizer(wx.HORIZONTAL)
      
        self.LeftPanel = wx.lib.scrolledpanel.ScrolledPanel(self, -1, style = wx.SUNKEN_BORDER, size = (500, 700))
        self.RightPanel = wx.lib.scrolledpanel.ScrolledPanel(self, -1, style = wx.SUNKEN_BORDER, size = (500, 700))
        
        self.LeftPanel.SetupScrolling(False, True)
        self.RightPanel.SetupScrolling(False, True)
        
        LeftBoxSizer = wx.BoxSizer(wx.VERTICAL)
        self.RightBoxSizer = wx.BoxSizer(wx.VERTICAL)
        
        BoxSizer.Add(self.LeftPanel, 1, wx.EXPAND|wx.ALL)
        BoxSizer.Add(self.RightPanel, 1, wx.EXPAND|wx.ALL)
        
        # Widgets for .card-file
        LeftBoxSizer.AddSpacer(10)
        
        CardBoxSizer = wx.BoxSizer(wx.HORIZONTAL)
        CardBoxSizer.AddSpacer(10)
        
        TextCard = wx.StaticText(self.LeftPanel, -1, 'Location of data:')
        CardBoxSizer.Add(TextCard)
        
        LeftBoxSizer.Add(CardBoxSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(10)
        
        CardPathSizer = wx.BoxSizer(wx.HORIZONTAL)
        CardPathSizer.AddSpacer(10)    
        
        self.DataPathTxt = wx.StaticText(self.LeftPanel, -1, 'N/A')
        CardPathSizer.Add(self.DataPathTxt)
        
        LeftBoxSizer.Add(CardPathSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(10)
        
        CardBtnBoxSizer = wx.BoxSizer(wx.HORIZONTAL)
        CardBtnBoxSizer.AddSpacer(10)
        
        BrowseDataBtn = wx.Button(self.LeftPanel, label = 'Browse')
        self.Bind(wx.EVT_BUTTON, self.BrowseDataFnc, BrowseDataBtn)
        CardBtnBoxSizer.Add(BrowseDataBtn, 1, wx.EXPAND)
        CardBtnBoxSizer.AddSpacer(10)
        
        LeftBoxSizer.Add(CardBtnBoxSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(10)
        
        LineCard = wx.StaticLine(self.LeftPanel)
        LeftBoxSizer.Add(LineCard, wx.EXPAND|wx.HORIZONTAL)
        
        self.DataPathStr = 'N/A'
        
        # Widgets for sample-file
        LeftBoxSizer.AddSpacer(10)
        
        SampleBoxSizer = wx.BoxSizer(wx.HORIZONTAL)
        SampleBoxSizer.AddSpacer(10)
        
        TextSamples = wx.StaticText(self.LeftPanel, -1, 'Location of sample info:')
        SampleBoxSizer.Add(TextSamples)
        
        LeftBoxSizer.Add(SampleBoxSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(10)
        
        SamplePathSizer = wx.BoxSizer(wx.HORIZONTAL)
        SamplePathSizer.AddSpacer(10)
        
        self.SamplesPathTxt = wx.StaticText(self.LeftPanel, -1, 'N/A')
        SamplePathSizer.Add(self.SamplesPathTxt)
        
        LeftBoxSizer.Add(SamplePathSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(10)
        
        SamplesBtnSizer = wx.BoxSizer(wx.HORIZONTAL)
        SamplesBtnSizer.AddSpacer(10)
        
        BrowseSamplesBtn = wx.Button(self.LeftPanel, label = 'Browse')
        self.Bind(wx.EVT_BUTTON, self.BrowseSamplesFnc, BrowseSamplesBtn)
        SamplesBtnSizer.Add(BrowseSamplesBtn, 1, wx.EXPAND)
        SamplesBtnSizer.AddSpacer(10)
        
        LeftBoxSizer.Add(SamplesBtnSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(10)
        
        LineSamples = wx.StaticLine(self.LeftPanel)
        LeftBoxSizer.Add(LineSamples, wx.EXPAND|wx.HORIZONTAL)
        
        self.SamplesPathStr = 'N/A'
        
        # Widgets for specifying q-range
        LeftBoxSizer.AddSpacer(10)
        
        MinimumQBoxSizer = wx.BoxSizer(wx.HORIZONTAL)
        MinimumQBoxSizer.AddSpacer(10)
        
        TextMinimumQ = wx.StaticText(self.LeftPanel, -1, 'Lower and upper limit in q:')
        MinimumQBoxSizer.Add(TextMinimumQ)
        
        LeftBoxSizer.Add(MinimumQBoxSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(2)  
        
        QBoxSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.QMinTxtCtrl = wx.TextCtrl(self.LeftPanel, -1, '', size = (80, -1), style = wx.TE_CENTRE)
        self.QMinTxtCtrl.SetValue('0')
        QBoxSizer.Add(self.QMinTxtCtrl)
        
        QBoxSizer.AddSpacer(20)
        
        self.QMaxTxtCtrl = wx.TextCtrl(self.LeftPanel, -1, '', size = (80, -1), style = wx.TE_CENTRE)
        self.QMaxTxtCtrl.SetValue('1')
        QBoxSizer.Add(self.QMaxTxtCtrl)
        
        LeftBoxSizer.Add(QBoxSizer, 0, wx.ALIGN_CENTER)
        LeftBoxSizer.AddSpacer(10)
        
        LineQ = wx.StaticLine(self.LeftPanel)
        LeftBoxSizer.Add(LineQ, wx.EXPAND|wx.HORIZONTAL)
        
        # Widgets used to control, which function is used for modelling
        LeftBoxSizer.AddSpacer(10)
        
        ModelTextBoxSizer = wx.BoxSizer(wx.HORIZONTAL)
        ModelTextBoxSizer.AddSpacer(10)        
        
        ModelText = wx.StaticText(self.LeftPanel, -1, 'Model:')
        ModelTextBoxSizer.Add(ModelText)
        
        (ModelFunctionNames, self.ModelFunctionPaths) = self.GetModels()
        (ModelFunctionNames, self.ModelFunctionPaths) = zip(*sorted(zip(ModelFunctionNames, self.ModelFunctionPaths)))
        
        LeftBoxSizer.Add(ModelTextBoxSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(10)
        
        self.ModelComboBox = wx.ComboBox(self.LeftPanel, -1, 'Pick a model function...' , size = (300, -1), choices = ModelFunctionNames, style = wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.ComboBoxSelectFnc)
        
        LeftBoxSizer.Add(self.ModelComboBox, 0, wx.ALIGN_CENTER)
        LeftBoxSizer.AddSpacer(10)
        
        LineModel = wx.StaticLine(self.LeftPanel)
        LeftBoxSizer.Add(LineModel, wx.EXPAND|wx.HORIZONTAL)
        
        # Widgets used to choose algorithm
        LeftBoxSizer.AddSpacer(10)
        
        AlgorithmTextBoxSizer = wx.BoxSizer(wx.HORIZONTAL)
        AlgorithmTextBoxSizer.AddSpacer(10)   
        
        AlgorithmText = wx.StaticText(self.LeftPanel, -1, 'Algorithm:')
        AlgorithmTextBoxSizer.Add(AlgorithmText)
        
        LeftBoxSizer.Add(AlgorithmTextBoxSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(1)
        
        AlgorithmBtn1Sizer = wx.BoxSizer(wx.HORIZONTAL)
        AlgorithmBtn1Sizer.AddSpacer(10) 
        
        self.RBComputeModel = wx.RadioButton(self.LeftPanel, -1, 'Compute model', style = wx.RB_GROUP)
        AlgorithmBtn1Sizer.Add(self.RBComputeModel)
        
        LeftBoxSizer.Add(AlgorithmBtn1Sizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(1)
        
        AlgorithmBtn2Sizer = wx.BoxSizer(wx.HORIZONTAL)
        AlgorithmBtn2Sizer.AddSpacer(10)
        
        self.RBSteepestDescent = wx.RadioButton(self.LeftPanel, -1, 'Levenberg-Marquardt')
        AlgorithmBtn2Sizer.Add(self.RBSteepestDescent)
        
        LeftBoxSizer.Add(AlgorithmBtn2Sizer, 0, wx.EXPAND|wx.HORIZONTAL)
        AlgorithmBtn3Sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.CBPrintCorrelationMatrix = wx.CheckBox(self.LeftPanel, -1, 'Print covariance matrix')
        AlgorithmBtn3Sizer.Add(self.CBPrintCorrelationMatrix)
        
        LeftBoxSizer.Add(AlgorithmBtn3Sizer, 0, wx.ALIGN_CENTER)
        LeftBoxSizer.AddSpacer(1)
        
        AlgorithmBtn4Sizer = wx.BoxSizer(wx.HORIZONTAL)
        AlgorithmBtn4Sizer.AddSpacer(10)
        
        self.RBGridsearch = wx.RadioButton(self.LeftPanel, -1, 'Gridsearch (LM) - # of cycles:')
        AlgorithmBtn4Sizer.Add(self.RBGridsearch)
        
        LeftBoxSizer.Add(AlgorithmBtn4Sizer, 0, wx.EXPAND|wx.HORIZONTAL)
        AlgorithmBtn5Sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.CyclesSpn = wx.SpinCtrl(self.LeftPanel, -1, "", size = (100, -1), min = 1, max = 10000, initial = 10)        
        AlgorithmBtn5Sizer.Add(self.CyclesSpn)
        
        LeftBoxSizer.Add(AlgorithmBtn5Sizer, 0, wx.ALIGN_CENTER)
        LeftBoxSizer.AddSpacer(1)
        
        AlgorithmBtn7Sizer = wx.BoxSizer(wx.HORIZONTAL)
        AlgorithmBtn7Sizer.AddSpacer(10)
        
        self.RBBFGS = wx.RadioButton(self.LeftPanel, -1, 'Broyden-Fletcher-Goldfarb-Shanno')
        AlgorithmBtn7Sizer.Add(self.RBBFGS)
        
        ## Currently unavailable ##
        self.RBBFGS.Disable()
        ## Currently unavailable ##
        
        LeftBoxSizer.Add(AlgorithmBtn7Sizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(1)
        
        AlgorithmBtn8Sizer = wx.BoxSizer(wx.HORIZONTAL)
        AlgorithmBtn8Sizer.AddSpacer(10)        
        
        self.RBBFGSGridsearch = wx.RadioButton(self.LeftPanel, -1, 'Gridsearch (BFGS) - # of cycles:')
        AlgorithmBtn8Sizer.Add(self.RBBFGSGridsearch)
        
        ## Currently unavailable ##
        self.RBBFGSGridsearch.Disable()
        ## Currently unavailable ##
        
        LeftBoxSizer.Add(AlgorithmBtn8Sizer, 0, wx.EXPAND|wx.HORIZONTAL)
        AlgorithmBtn9Sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.CyclesBFGSSpn = wx.SpinCtrl(self.LeftPanel, -1, "", size = (100, -1), min = 1, max = 10000, initial = 10)
        AlgorithmBtn9Sizer.Add(self.CyclesBFGSSpn)
        
        ## Currently unavailable ##
        self.CyclesBFGSSpn.Disable()
        ## Currently unavailable ##
        
        LeftBoxSizer.Add(AlgorithmBtn9Sizer, 0, wx.ALIGN_CENTER)
        LeftBoxSizer.AddSpacer(1)
        
        AlgorithmBtn10Sizer = wx.BoxSizer(wx.HORIZONTAL)
        AlgorithmBtn10Sizer.AddSpacer(10)
        
        self.RBSwarm = wx.RadioButton(self.LeftPanel, -1, 'Swarm - # of agents - # of flights:')
        AlgorithmBtn10Sizer.Add(self.RBSwarm)
        
        LeftBoxSizer.Add(AlgorithmBtn10Sizer, 0, wx.EXPAND|wx.HORIZONTAL)
        SwarmSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.SwarmAgentsSpn = wx.SpinCtrl(self.LeftPanel, -1, "", size = (100, -1), min = 1, max = 10000, initial = 32)
        SwarmSizer.Add(self.SwarmAgentsSpn)
        
        SwarmSizer.AddSpacer(20)
        
        self.SwarmFlightsSpn = wx.SpinCtrl(self.LeftPanel, -1, "", size = (100, -1), min = 1, max = 10000, initial = 3)
        SwarmSizer.Add(self.SwarmFlightsSpn)
        
        LeftBoxSizer.Add(SwarmSizer, 0, wx.ALIGN_CENTER)
        LeftBoxSizer.AddSpacer(1)
        
        AlgorithmBtn11Sizer = wx.BoxSizer(wx.HORIZONTAL)
        AlgorithmBtn11Sizer.AddSpacer(10)
        
        self.RBGenetic = wx.RadioButton(self.LeftPanel, -1, 'Genetic - # of individuals - # of generations:')
        AlgorithmBtn11Sizer.Add(self.RBGenetic)
        
        LeftBoxSizer.Add(AlgorithmBtn11Sizer, 0, wx.EXPAND|wx.HORIZONTAL)
        GeneticSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.GeneticIndivudialsSpn = wx.SpinCtrl(self.LeftPanel, -1, "", size = (100, -1), min = 1, max = 10000, initial = 100)
        GeneticSizer.Add(self.GeneticIndivudialsSpn)
        
        GeneticSizer.AddSpacer(20)
        
        self.GeneticGenerationsSpn = wx.SpinCtrl(self.LeftPanel, -1, "", size = (100, -1), min = 1, max = 10000, initial = 10)
        GeneticSizer.Add(self.GeneticGenerationsSpn)
        
        LeftBoxSizer.Add(GeneticSizer, 0, wx.ALIGN_CENTER)
        LeftBoxSizer.AddSpacer(10)
        
        AlgorithmBtn6Sizer = wx.BoxSizer(wx.HORIZONTAL)
        AlgorithmBtn6Sizer.AddSpacer(10)
        
        self.RBProfileLikelihood = wx.RadioButton(self.LeftPanel, -1, 'Profile Likelihood - Fractile - # of cycles:')
        AlgorithmBtn6Sizer.Add(self.RBProfileLikelihood)
        
        LeftBoxSizer.Add(AlgorithmBtn6Sizer, 0, wx.EXPAND|wx.HORIZONTAL)
        ProfileLikelihoodSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.Chi2Ctrl = wx.TextCtrl(self.LeftPanel, -1, '', (100, -1), style = wx.TE_CENTRE)
        self.Chi2Ctrl.SetValue('3.84')
        ProfileLikelihoodSizer.Add(self.Chi2Ctrl)
        
        ProfileLikelihoodSizer.AddSpacer(20)
        
        self.CyclesPLHSpn = wx.SpinCtrl(self.LeftPanel, -1, "", size = (100, -1), min = 1, max = 10000, initial = 3)
        ProfileLikelihoodSizer.Add(self.CyclesPLHSpn)
        
        LeftBoxSizer.Add(ProfileLikelihoodSizer, 0, wx.ALIGN_CENTER)
        LeftBoxSizer.AddSpacer(1)
        
        AlgorithmBtn11Sizer = wx.BoxSizer(wx.HORIZONTAL)
        AlgorithmBtn11Sizer.AddSpacer(10)
        
        self.RBProfileLikelihoodSingle = wx.RadioButton(self.LeftPanel, -1, 'PLH (Single parameter) - Fractile - # of cycles - ID of parameter:')
        AlgorithmBtn11Sizer.Add(self.RBProfileLikelihoodSingle)
        
        LeftBoxSizer.Add(AlgorithmBtn11Sizer, 0, wx.EXPAND|wx.HORIZONTAL)
        ProfileLikelihoodSingleSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.Chi2CtrlSingle = wx.TextCtrl(self.LeftPanel, -1, '', (100, -1), style = wx.TE_CENTRE)
        self.Chi2CtrlSingle.SetValue('3.84')
        ProfileLikelihoodSingleSizer.Add(self.Chi2CtrlSingle)
        
        ProfileLikelihoodSingleSizer.AddSpacer(20)
        
        self.CyclesPLHSpnSingle = wx.SpinCtrl(self.LeftPanel, -1, "", size = (100, -1), min = 1, max = 10000, initial = 3)
        ProfileLikelihoodSingleSizer.Add(self.CyclesPLHSpnSingle)
        
        ProfileLikelihoodSingleSizer.AddSpacer(20)
        
        self.IDSinglePLHSpn = wx.SpinCtrl(self.LeftPanel, -1, "", size = (100, -1), min = 0, max = 10000, initial = 0)
        ProfileLikelihoodSingleSizer.Add(self.IDSinglePLHSpn)
        
        LeftBoxSizer.Add(ProfileLikelihoodSingleSizer, 0, wx.ALIGN_CENTER)
        LeftBoxSizer.AddSpacer(1)
        
        AlgorithmLine = wx.StaticLine(self.LeftPanel)
        LeftBoxSizer.Add(AlgorithmLine, wx.EXPAND|wx.HORIZONTAL)
        
        # Widgets used to control resolution function options
        LeftBoxSizer.AddSpacer(10)
        
        ResolutionTextBoxSizer = wx.BoxSizer(wx.HORIZONTAL)
        ResolutionTextBoxSizer.AddSpacer(10)  
        
        ResolutionText = wx.StaticText(self.LeftPanel, -1, 'Resolution effects:')
        ResolutionTextBoxSizer.Add(ResolutionText)
        
        LeftBoxSizer.Add(ResolutionTextBoxSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(1)
        
        Reso1Sizer = wx.BoxSizer(wx.HORIZONTAL)
        Reso1Sizer.AddSpacer(10)
        
        self.RBNoResolution = wx.RadioButton(self.LeftPanel, -1, 'Do not include resolution effects', style = wx.RB_GROUP)
        Reso1Sizer.Add(self.RBNoResolution)
        
        LeftBoxSizer.Add(Reso1Sizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(1)
        
        Reso2Sizer = wx.BoxSizer(wx.HORIZONTAL)
        Reso2Sizer.AddSpacer(10)
        
        self.RBResolution = wx.RadioButton(self.LeftPanel, -1, 'Include resolution effects - # of steps:')  
        Reso2Sizer.Add(self.RBResolution)
        
        LeftBoxSizer.Add(Reso2Sizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(1)
       
        Reso3Sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.ResolutionSpn = wx.SpinCtrl(self.LeftPanel, -1, "", size = (100, -1), min = 5, max = 25, initial = 15)
        Reso3Sizer.Add(self.ResolutionSpn)
        
        LeftBoxSizer.Add(Reso3Sizer, 0, wx.ALIGN_CENTER)
        LeftBoxSizer.AddSpacer(10)
        
        ResolutionBoxSizer = wx.BoxSizer(wx.HORIZONTAL)
        ResolutionBoxSizer.AddSpacer(10)
        
        ResolutionLocationText = wx.StaticText(self.LeftPanel, -1, 'Location of .res-file:')
        ResolutionBoxSizer.Add(ResolutionLocationText)
        
        LeftBoxSizer.Add(ResolutionBoxSizer, 0, wx.EXPAND|wx.HORIZONTAL)   
        LeftBoxSizer.AddSpacer(10)
        
        ResolutionPathBoxSizer = wx.BoxSizer(wx.HORIZONTAL)
        ResolutionPathBoxSizer.AddSpacer(10)
        
        self.ResolutionPathTxt = wx.StaticText(self.LeftPanel, -1, 'N/A')
        ResolutionPathBoxSizer.Add(self.ResolutionPathTxt)
        
        LeftBoxSizer.Add(ResolutionPathBoxSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(10)
        
        ResoBtnBoxSizer = wx.BoxSizer(wx.HORIZONTAL)
        ResoBtnBoxSizer.AddSpacer(10)
        
        BrowseResolutionBtn = wx.Button(self.LeftPanel, label = 'Browse')
        self.Bind(wx.EVT_BUTTON, self.BrowseResolutionFnc, BrowseResolutionBtn)
        ResoBtnBoxSizer.Add(BrowseResolutionBtn, 1, wx.EXPAND)
        ResoBtnBoxSizer.AddSpacer(10)
        
        LeftBoxSizer.Add(ResoBtnBoxSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(10)
        
        ResolutionLine = wx.StaticLine(self.LeftPanel, -1)
        LeftBoxSizer.Add(ResolutionLine, wx.EXPAND|wx.HORIZONTAL)
        
        # Widgets used to input .pdb-file if needed
        LeftBoxSizer.AddSpacer(10)
        
        PDBTextBoxSizer = wx.BoxSizer(wx.HORIZONTAL)
        PDBTextBoxSizer.AddSpacer(10)
        
        PDBText = wx.StaticText(self.LeftPanel, -1, 'Location of .pdb-file:')
        PDBTextBoxSizer.Add(PDBText)
        
        LeftBoxSizer.Add(PDBTextBoxSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(10)
        
        PDBPathSizer = wx.BoxSizer(wx.HORIZONTAL)
        PDBPathSizer.AddSpacer(10)
        
        self.PDBPathTxt = wx.StaticText(self.LeftPanel, -1, 'N/A')
        PDBPathSizer.Add(self.PDBPathTxt)
        
        LeftBoxSizer.Add(PDBPathSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(10)
        
        PDBBtnSizer = wx.BoxSizer(wx.HORIZONTAL)
        PDBBtnSizer.AddSpacer(10)
        
        BrowsePDBBtn = wx.Button(self.LeftPanel, label = 'Browse')
        self.Bind(wx.EVT_BUTTON, self.BrowsePDBFnc, BrowsePDBBtn)
        self.PDBPathStr = 'N/A'
        PDBBtnSizer.Add(BrowsePDBBtn, 1, wx.EXPAND)
        PDBBtnSizer.AddSpacer(10)
        
        LeftBoxSizer.Add(PDBBtnSizer, 0, wx.EXPAND|wx.HORIZONTAL)
        LeftBoxSizer.AddSpacer(10)
        
        PDBLine = wx.StaticLine(self.LeftPanel)
        LeftBoxSizer.Add(PDBLine, wx.EXPAND|wx.HORIZONTAL)
        
        # Widgets used on the right side of the panel
        self.RightBoxSizer.AddSpacer(10)
        
        ParBtnSpacer = wx.BoxSizer(wx.HORIZONTAL)
        ParBtnSpacer.AddSpacer(10)
        
        self.BrowseParametersBtn = wx.Button(self.RightPanel, label = 'Load parameter file')
        self.Bind(wx.EVT_BUTTON, self.BrowseParametersFnc, self.BrowseParametersBtn)
        self.ParametersCreated = False
        self.BrowseParametersBtn.Disable()
        
        ParBtnSpacer.Add(self.BrowseParametersBtn, 1, wx.EXPAND)
        ParBtnSpacer.AddSpacer(10)
        
        self.RightBoxSizer.Add(ParBtnSpacer, 0, wx.EXPAND|wx.HORIZONTAL)
        self.RightBoxSizer.AddSpacer(5)
        
        CategorySizer = wx.BoxSizer(wx.HORIZONTAL)
        CategorySizer.AddSpacer(10)
        
        CategorySizer.Add(wx.StaticText(self.RightPanel, -1, 'Minimum:', size = (70, -1)))
        CategorySizer.AddSpacer(20)
        CategorySizer.Add(wx.StaticText(self.RightPanel, -1, 'Value:', size = (100, -1)))
        CategorySizer.AddSpacer(20)
        CategorySizer.Add(wx.StaticText(self.RightPanel, -1, 'Maximum:', size = (70, -1)))
        CategorySizer.AddSpacer(20)
        CategorySizer.Add(wx.StaticText(self.RightPanel, -1, 'Fit?:', size = (30, -1)))
        CategorySizer.AddSpacer(20)
        CategorySizer.Add(wx.StaticText(self.RightPanel, -1, 'Name:', size = (100, -1)))
        
        self.RightBoxSizer.Add(CategorySizer, 0, wx.EXPAND|wx.HORIZONTAL)
        self.RightBoxSizer.AddSpacer(2)
        
        # Main Control Buttons
        LeftBoxSizer.AddSpacer(10) 
        
        ButtonBoxSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.ExecuteBtn = wx.Button(self.LeftPanel, label = 'Execute')
        self.Bind(wx.EVT_BUTTON, self.ExecuteFnc, self.ExecuteBtn)
        self.ExecuteBtn.Disable()
        ButtonBoxSizer.Add(self.ExecuteBtn, 0, wx.EXPAND|wx.HORIZONTAL)
        
        ButtonBoxSizer.AddSpacer(10)
        
        self.UndoBtn = wx.Button(self.LeftPanel, label = 'Undo last')
        self.Bind(wx.EVT_BUTTON, self.UndoFnc, self.UndoBtn)
        self.UndoBtn.Disable()
        ButtonBoxSizer.Add(self.UndoBtn, 0, wx.EXPAND|wx.HORIZONTAL)
        
        ButtonBoxSizer.AddSpacer(10)
        
        self.PlotBtn = wx.Button(self.LeftPanel, label = 'Plot')
        self.Bind(wx.EVT_BUTTON, self.PlotFnc, self.PlotBtn)
        self.PlotBtn.Disable()
        ButtonBoxSizer.Add(self.PlotBtn, 0, wx.EXPAND|wx.HORIZONTAL)
        
        ButtonBoxSizer.AddSpacer(10)
        
        QuitBtn = wx.Button(self.LeftPanel, label = 'Quit')
        self.Bind(wx.EVT_BUTTON, self.CloseWindowFnc, QuitBtn)
        self.Bind(wx.EVT_CLOSE, self.CloseWindowFnc)
        ButtonBoxSizer.Add(QuitBtn, 0, wx.EXPAND|wx.HORIZONTAL)
        
        LeftBoxSizer.Add(ButtonBoxSizer, 0, wx.ALIGN_CENTER)
        
        # Conclusion of init-function
        self.LeftPanel.SetSizer(LeftBoxSizer)
        self.RightPanel.SetSizer(self.RightBoxSizer)
        self.SetSizer(BoxSizer)
        self.SetAutoLayout(True)
        self.LeftPanel.Layout()
        self.RightPanel.Layout()
        self.Layout()
        
        # Check for existence of compiler
        (CompilerPath, Compilerflags) = self.SetupCompiler()
        
        try:
            CompilerCheck = subprocess.call([CompilerPath, "-v"])
        except:
            CompilerCheck = -1
        
        if CompilerCheck == 0:
            print ("")
            print ("")
            print ("")
            print ("************************")
            print ("* Welcome to WillItFit *")
            print ("************************")
            print ("")
            print ("Found all necessary Python-libraries.")
            print ("Found compiler: " + CompilerPath + ".")
        else:
            wx.MessageBox("Unable to locate the C-compiler specified in CompilerConfig.cfg. This compiler must be installed and correctly set up in order for WillItFit to work.", \
                          "No compiler found...!", wx.OK | wx.ICON_INFORMATION)

## Define function to be executed
    def ExecuteFnc(self, event):

        self.ClearOldData()
        os.mkdir('./.data/')
        
        # Get paramters from different widgets
        ChiSquareFractile = 0.0
        MinQStr = self.QMinTxtCtrl.GetValue()
        MaxQStr = self.QMaxTxtCtrl.GetValue()
        DataStr = self.DataPathStr
        SamplesStr = self.SamplesPathStr
        PDBStr = self.PDBPathStr
        
        if self.RBComputeModel.GetValue():
            FittingRoutine = 0
            FittingRoutineArgument2 = 0
            FittingRoutineArgument3 = 0
        elif self.RBSteepestDescent.GetValue():
            FittingRoutine = 1
            FittingRoutineArgument2 = 0
            FittingRoutineArgument3 = 0
        elif self.RBGridsearch.GetValue():
            FittingRoutine = 2
            FittingRoutineArgument2 = self.CyclesSpn.GetValue()
            FittingRoutineArgument3 = 0
        elif self.RBBFGS.GetValue():
            FittingRoutine = 3
            FittingRoutineArgument2 = 0
            FittingRoutineArgument3 = 0
        elif self.RBBFGSGridsearch.GetValue():
            FittingRoutine = 4
            FittingRoutineArgument2 = self.CyclesBFGSSpn.GetValue()
            FittingRoutineArgument3 = 0
        elif self.RBSwarm.GetValue():
            FittingRoutine = 5
            FittingRoutineArgument2 = self.SwarmFlightsSpn.GetValue()
            FittingRoutineArgument3 = self.SwarmAgentsSpn.GetValue()
        elif self.RBGenetic.GetValue():
            FittingRoutine = 6
            FittingRoutineArgument2 = self.GeneticGenerationsSpn.GetValue()
            FittingRoutineArgument3 = self.GeneticIndivudialsSpn.GetValue()
        elif self.RBProfileLikelihood.GetValue():
            FittingRoutine = 7
            ChiSquareFractile = self.Chi2Ctrl.GetValue()
            FittingRoutineArgument2 = 0
            FittingRoutineArgument3 = self.CyclesPLHSpn.GetValue()
        else:
            FittingRoutine = 8
            ChiSquareFractile = self.Chi2CtrlSingle.GetValue()
            FittingRoutineArgument2 = self.IDSinglePLHSpn.GetValue()
            FittingRoutineArgument3 = self.CyclesPLHSpnSingle.GetValue()
        
        if self.RBNoResolution.GetValue():
            NumberOfSteps = 0
            ResolutionFileName = 'N/A'
        else:
            NumberOfSteps = self.ResolutionSpn.GetValue()
            ResolutionFileName = self.ResolutionPathStr
        
        if self.CBPrintCorrelationMatrix.GetValue():
            PrintCorrelationInt = 1
        else:
            PrintCorrelationInt = 0
        
        # Build dummy parameter file
        file = open('ParametersBefore.mcp','w')
        
        i = 0
        OnOffToPrint = []
        NoFreeParameters = True
        
        for line in self.ParameterValues:
            if self.Checkboxes[i].GetValue():
                OnOffToPrint.append('0')
                NoFreeParameters = False
            else:
                OnOffToPrint.append('1')
            
            LineToWrite = self.ParameterMinTxt[i].GetValue() + '   ' + \
                          self.ParameterValuesTxt[i].GetValue() + '   ' + \
                          self.ParameterMaxTxt[i].GetValue() + '   ' + \
                          OnOffToPrint[i] + '   ' + \
                          self.ParameterNames[i] + '\n'
            
            file.write(LineToWrite)
            
            i = i + 1
        
        file.close()
        
        if NoFreeParameters and (FittingRoutine != 0):
            wx.MessageBox('This fitting routine needs free parameters...', "It won't fit...", wx.OK | wx.ICON_INFORMATION)
            return
        
        ParametersStr = 'ParametersBefore.mcp'
        
        try:
            os.remove('ParametersAfter.mcp')
        except:
            pass
        
        # Prepare algorithm
        Cardfilename = os.path.basename(DataStr)
        
        try:
            self.ResultsDirectory = './' + Cardfilename + '-results/'
            os.mkdir(self.ResultsDirectory)
        except:
            pass
        
        if FittingRoutine == 7 or FittingRoutine == 8:
            try:
                os.mkdir('ProfileLikelihood')
            except:
                pass
        
        # Compile C-module and run
        SourceFileLocation = os.path.join(self.InitialDirectoryStr, "WillItFit.c")
        ExecutableFileLocation = os.path.join(self.InitialDirectoryStr, "WillItFit.com")
        
        (CompilerPath, Compilerflags) = self.SetupCompiler()
        
        if "-fopenmp" in Compilerflags:
            print ("Using OpenMP...")            
            File = open(os.path.join(self.InitialDirectoryStr, "Auxillary", "Parallelisation.h"), "w+")
            File.write("#include <omp.h>")
            File.close()
        else:
            File = open(os.path.join(self.InitialDirectoryStr, "Auxillary", "Parallelisation.h"), "w+")
            File.write(" ")
            File.close()
        
        if "UserDefined.h" in os.listdir(os.path.join(self.InitialDirectoryStr, 'Models', self.ModelPath)):
            print ("Custom-made UserDefined.h located...")
            File = open(os.path.join(self.InitialDirectoryStr, "Auxillary", "IncludeUserDefined.h"), "w+")
            File.write('#include "' + os.path.join(self.InitialDirectoryStr, 'Models', self.ModelPath, 'UserDefined.h'))
            File.close()
        else:
            File = open(os.path.join(self.InitialDirectoryStr, "Auxillary", "IncludeUserDefined.h"), "w+")
            File.write('#include "UserDefined.h"')
            File.close()
        
        if "WillItFit.com" not in os.listdir(self.InitialDirectoryStr):
            print ("Compiling (no compiled C-module found...)")          
            CompilerArguments = [CompilerPath, SourceFileLocation, "-o", ExecutableFileLocation] + Compilerflags
            print (CompilerArguments)
            print ("\n")
            Compilationresult = subprocess.call(CompilerArguments)
        elif os.path.getmtime(ExecutableFileLocation) < self.CheckNewestFile():
            os.remove(ExecutableFileLocation)
            print ("Recompiling (the source code has been altered, or a new model has been selected...)")
            CompilerArguments = [CompilerPath, SourceFileLocation, "-o", ExecutableFileLocation] + Compilerflags
            print (CompilerArguments)
            print ("\n")
            Compilationresult = subprocess.call(CompilerArguments)
        else:
            Compilationresult = 0        
        
        if Compilationresult != 0:
            wx.MessageBox("The compilation of WillItFit failed. \nThe system tried to use " + CompilerPath + " with the following arguments: \n" + \
                          " ".join([str(x) for x in Compilerflags]) + ". \nAre all of these available? Alternatively, there might be an error in the syntax.", \
                          "It cannot fit...?", wx.OK | wx.ICON_INFORMATION)
            return
        
        if "WillItFit.com" not in os.listdir(self.InitialDirectoryStr):
            wx.MessageBox("Unable to find the compiled WillItFit.com-executable. Was the program compiled correctly?", "It cannot fit...?", wx.OK | wx.ICON_INFORMATION)
            return
        
        ProcessToCall = []
        ProcessToCall.append(ExecutableFileLocation)
        ProcessToCall.append('-c=%s' % DataStr)
        ProcessToCall.append('-s=%s' % SamplesStr)
        ProcessToCall.append('-p=%s' % ParametersStr)
        ProcessToCall.append('-n=%s' % MinQStr)
        ProcessToCall.append('-x=%s' % MaxQStr)
        ProcessToCall.append('-r=%d' % FittingRoutine)
        ProcessToCall.append('-y=%d' % FittingRoutineArgument2)
        ProcessToCall.append('-t=%d' % NumberOfSteps)
        ProcessToCall.append('-e=%s' % ResolutionFileName)
        ProcessToCall.append('-o=%d' % PrintCorrelationInt)
        ProcessToCall.append('-d=%s' % PDBStr)
        ProcessToCall.append('-h=%s' % ChiSquareFractile)
        ProcessToCall.append('-a=%d' % FittingRoutineArgument3)
        print (ProcessToCall)
        print ("\n")

        Process = subprocess.Popen(ProcessToCall)
        BusyDialog = WaitingDialog("WillItFit is running...", "Will it fit...?", Process)
        if BusyDialog.ShowModal() == wx.ID_CANCEL:
            Process.kill()
            DisplayedMessage = 'Algorithm aborted prematurely...'
            wx.MessageBox(DisplayedMessage, "It didn't fit...?", wx.OK | wx.ICON_INFORMATION)
            return
        
        # Enable buttons
        self.UndoBtn.Enable()
        self.PlotBtn.Enable()
        
        # Display concluding message and new parameters
        try:
            self.PrintParameters('ParametersAfter.mcp')
            self.PrintResultsToCardfileDirectory()
            
            Chisquare = GetChisquare()
            DisplayedMessage = 'Algorithm concluded correctly and returned a chisquare of ' + Chisquare + '.'
            wx.MessageBox(DisplayedMessage, 'Info', wx.OK | wx.ICON_INFORMATION)      
        except:
            ReturnMessage = GetReturnMessage()
            DisplayedMessage = 'Algorithm concluded incorrectly with the following message: ' + ReturnMessage
            wx.MessageBox(DisplayMessage, 'Info', wx.OK | wx.ICON_INFORMATION)
## Define function to be executed
    def UndoFnc(self, event):
        self.PrintParameters('ParametersBefore.mcp')
        self.UndoBtn.Disable()

## Define plotting function
    def PlotFnc(self, event):
        pylab.clf()
        
        DataArr = []
        FitArr = [] 
        
        for Filename in os.listdir('./.data/'):
            if os.path.splitext(Filename)[1] == '.dat':
                ShortName = os.path.basename(Filename)
                
                if ShortName[0:3] == 'fit':
                    FitArr.append(os.path.join('.data/', ShortName))
                    
                if ShortName[0:3] == 'dat':
                    if ShortName[0:4] != 'data':
                        DataArr.append(os.path.join('.data/', ShortName))
                                       
        NumberOfSets = len(DataArr)
        
        DataArr.sort()
        FitArr.sort()
        
        xData = []
        yData = []
        wData = []
        
        xFit = []
        yFit = []
        
        ContrastData = []
        ContrastFit = []
        
        for i in range(NumberOfSets):
            xData.append([])
            yData.append([])
            wData.append([])
            
            xFit.append([])
            yFit.append([])
        
        for i in range(NumberOfSets):
            xDataArr = GetData(DataArr[i])[0]
            yDataArr = GetData(DataArr[i])[1]
            wDataArr = GetData(DataArr[i])[2]
            
            ContrastData.append(GetData(DataArr[i])[3])
            NumberOfPoints = len(xDataArr)
            
            for j in range(NumberOfPoints):
                if yDataArr[j] > 0:
                    xData[i].append(xDataArr[j])
                    yData[i].append(yDataArr[j])
                    wData[i].append(wDataArr[j])                   
            
            xFitArr = GetData(FitArr[i])[0]
            yFitArr = GetData(FitArr[i])[1]
            
            ContrastFit.append(GetData(DataArr[i])[3])
            NumberOfPoints = len(xFitArr)
            
            for j in range(NumberOfPoints):
                if yFitArr[j] > 0:
                    xFit[i].append(xFitArr[j])
                    yFit[i].append(yFitArr[j])
        
        Figure = pylab.figure(1) 
        Subplot = pylab.subplot(111)

        LowestYValue = 100000000000000
        
        for i in range(NumberOfSets):
            Subplot.set_xscale('log', nonposx = 'clip')
            Subplot.set_yscale('log', nonposy = 'clip')
            Subplot.errorbar(xData[i], yData[i], yerr = wData[i], fmt = FormatFnc(ContrastData[i]), label = 'Data' + str(i))
            
            for y in yData[i]:
                if y < LowestYValue:
                    LowestYValue = y

        for i in range(NumberOfSets):
            Subplot.set_xscale('log', nonposx = 'clip')
            Subplot.set_yscale('log', nonposy = 'clip')
            Subplot.plot(xFit[i], yFit[i], color = 'black')
        
        Subplot.set_ylim(ymin = LowestYValue / 10.0)
        
        Subplot.set_xlabel('q / (1/A)')
        Subplot.set_ylabel('I / (1/cm)')
        
        pylab.show()

## Define function bound to combobox
    def EnableExecute(self):
        if self.DataPathStr != 'N/A':
            if self.SamplesPathStr != 'N/A':
                if self.ParametersCreated:
                    self.ExecuteBtn.Enable()
    
    def ComboBoxSelectFnc(self, event):
        self.ModelFunctionInt = event.GetSelection()
        self.ModelPath = self.ModelFunctionPaths[self.ModelFunctionInt]
        
        ParameterFileFound = False
        
        HeaderFile = open(os.path.join(self.InitialDirectoryStr, 'Auxillary', 'ModelLocation.h'), "w+")
        HeaderFile.write('#include "' + os.path.join(self.InitialDirectoryStr, 'Models', self.ModelPath, 'ModelInfo.h') + '"')
        HeaderFile.close()
        
        for File in os.listdir(os.path.join(self.InitialDirectoryStr, 'Models', self.ModelPath)):
            if os.path.splitext(File)[1] == '.par':
                ParameterFileName = File
                ParameterFileFound = True
        
        if ParameterFileFound == False:
            wx.MessageBox("Unable to locate .par-file - is it present in the model directory?", "It cannot fit...?", wx.OK | wx.ICON_INFORMATION)
            return
        
        if self.ParametersCreated:
            self.RemoveParameters()
        
        self.ParametersCreated = True        
        PathToParameterFile = os.path.join(os.path.join(self.InitialDirectoryStr, 'Models', self.ModelPath, ParameterFileName))
        
        self.CreateParameters(PathToParameterFile)
        self.PrintParameters(PathToParameterFile)
        
        self.EnableExecute()
        self.BrowseParametersBtn.Enable()
        
    def CheckNewestFile(self):
        NewestTime = 0
        
        for Path, Subfolder, Files in os.walk(self.InitialDirectoryStr):
            for Filename in Files:
                if Filename not in ("Parallelisation.h", "WillItFit.com", "IncludeUserDefined.h"):
                    if "Examples" not in Path:
                        PathToFile = os.path.join(Path, Filename)
                        TimeOfModification = float(os.path.getmtime(PathToFile))
                        
                        if NewestTime < TimeOfModification:
                            NewestTime = TimeOfModification
        
        return NewestTime

## Define function used to search for viable models
    def GetModels(self):
        ListOfModelPaths = os.listdir(os.path.join(self.InitialDirectoryStr, 'Models'))
        ListOfModelNames = []
        
        for Model in ListOfModelPaths:
            if Model == '.DS_Store':
                Modelname = 'ZZ NO model! .DS_Store (ZZ to be last in list)'
                ListOfModelNames.append(Modelname.strip())
            if Model != '.DS_Store':
                Infofile = open(os.path.join(self.InitialDirectoryStr, 'Models', Model, 'ModelInfo.h'), 'rU')
                ModelInfoTxt = Infofile.readlines()
                Infofile.close()
            
                for Line in ModelInfoTxt:
                    try:
                        Modelname = re.split(r' * Model name: ', Line)[1]
                    except:
                        Modelname = False
                
                    if Modelname:
                        ListOfModelNames.append(Modelname.strip())
        
        return ListOfModelNames, ListOfModelPaths
    
## Define function for filebrowsing for data
    def BrowseDataFnc(self, event):
        FileDialogWindow = wx.FileDialog(None, 'Please select .card-file...', os.getcwd(), defaultFile = '', \
                                         wildcard = 'Card files (.card) |*.card|' 'All files |*.', style = wx.FD_OPEN)
        
        if FileDialogWindow.ShowModal() == wx.ID_OK:
            self.DataPathStr = FileDialogWindow.GetPath()
            DataPathDisplayStr = str(self.DataPathStr)
            
            while len(DataPathDisplayStr) > 40:
                DataPathDisplayStr = DataPathDisplayStr[1:]
            
            if len(self.DataPathStr) > 40:
                DataPathDisplayStr = '...' + DataPathDisplayStr
            
            self.DataPathTxt.SetLabel(DataPathDisplayStr)
            
            DirectoryStr = os.path.dirname(self.DataPathStr)            
            os.chdir(DirectoryStr)
            
            self.UndoBtn.Disable()
            self.PlotBtn.Disable()
            self.EnableExecute()
        
        FileDialogWindow.Destroy()
        
        # Clean up
        self.ClearOldParameters()

## Define function for filebrowsing for samples
    def BrowseSamplesFnc(self, event):
        FileDialogWindow = wx.FileDialog(None, 'Please select samples file...', os.getcwd(), defaultFile = '', \
                                         wildcard = 'Sample files (.dat) |*.dat|' 'All files |*.', style = wx.FD_OPEN)
        
        if FileDialogWindow.ShowModal() == wx.ID_OK:
            self.SamplesPathStr = FileDialogWindow.GetPath()
            SamplesPathDisplayStr = str(self.SamplesPathStr)
            
            while len(SamplesPathDisplayStr) > 40:
                SamplesPathDisplayStr = SamplesPathDisplayStr[1:]
            
            if len(self.SamplesPathStr) > 40:
                SamplesPathDisplayStr = '...' + SamplesPathDisplayStr
            
            self.SamplesPathTxt.SetLabel(SamplesPathDisplayStr)
            
            self.UndoBtn.Disable()
            self.PlotBtn.Disable()
            self.EnableExecute()
        
        FileDialogWindow.Destroy()

## Define function for filebrowsing for parameters
    def BrowseParametersFnc(self, event):
        FileDialogWindow = wx.FileDialog(None, 'Please select parameter file...', os.getcwd(), defaultFile = '', \
                                         wildcard = 'Parameter files (.par) |*.par|' 'All files |*.', style = wx.FD_OPEN)
        
        if FileDialogWindow.ShowModal() == wx.ID_OK:
            self.ParametersPathStr = FileDialogWindow.GetPath()
            
            if self.ParametersCreated:
                self.RemoveParameters()
            
            self.ParametersCreated = True
            self.CreateParameters(self.ParametersPathStr)
            self.PrintParameters(self.ParametersPathStr)
            
            self.UndoBtn.Disable()
            self.PlotBtn.Disable()
            self.EnableExecute()
        
        FileDialogWindow.Destroy()

## Define function for filebrowsing for samples
    def BrowseResolutionFnc(self, event):
        FileDialogWindow = wx.FileDialog(None, 'Please select resolution info file...', os.getcwd(), defaultFile = '', \
                                         wildcard = 'Resolution files (.res) |*.res|' 'All files |*.', style = wx.FD_OPEN)
        
        if FileDialogWindow.ShowModal() == wx.ID_OK:
            self.ResolutionPathStr = FileDialogWindow.GetPath()
            ResolutionPathDisplayStr = str(self.ResolutionPathStr)
            
            while len(ResolutionPathDisplayStr) > 40:
                ResolutionPathDisplayStr = ResolutionPathDisplayStr[1:]
            
            if len(self.ResolutionPathStr) > 40:
                ResolutionPathDisplayStr = '...' + ResolutionPathDisplayStr
            
            self.ResolutionPathTxt.SetLabel(ResolutionPathDisplayStr)
            
            self.UndoBtn.Disable()
            self.PlotBtn.Disable()
        
        FileDialogWindow.Destroy()
        
## Define function for filebrowsing for samples
    def BrowsePDBFnc(self, event):
        FileDialogWindow = wx.FileDialog(None, 'Please select .pdb-file...', os.getcwd(), defaultFile = '', \
                                         wildcard = 'PDB-files (.pdb) |*.pdb|' 'All files |*.', style = wx.FD_OPEN)
        
        if FileDialogWindow.ShowModal() == wx.ID_OK:
            self.PDBPathStr = FileDialogWindow.GetPath()
            PDBPathDisplayStr = str(self.PDBPathStr)
            
            while len(PDBPathDisplayStr) > 40:
                PDBPathDisplayStr = PDBPathDisplayStr[1:]
            
            if len(self.PDBPathStr) > 40:
                PDBPathDisplayStr = '...' + PDBPathDisplayStr
            
            self.PDBPathTxt.SetLabel(PDBPathDisplayStr)
            
            self.UndoBtn.Disable()
            self.PlotBtn.Disable()
        
        FileDialogWindow.Destroy()
        
## Compilersetup
    def SetupCompiler(self):
        CompilerInfoFile = open(os.path.join(self.InitialDirectoryStr, 'CompilerInfo.cfg'), 'rU')
        CompilerInfo = CompilerInfoFile.readlines()
        CompilerInfoFile.close()
        
        for Line in CompilerInfo:
            try:
                CompilerPath = re.split(r'CompilerPath  = ', Line)[1].strip()
            except:
                pass
            
            try:
                Compilerflags = re.split(r'Compilerflags = ', Line)[1].strip()
            except:
                pass
        
        Compilerflags = re.split(" ", Compilerflags)
        
        return CompilerPath, Compilerflags

## Print the obtained parameters, names and values
    def CreateParameters(self, filename):
        file = open(filename, 'rU')
        self.ParametersArr = file.readlines()
        file.close()
        
        NumberOfParameters = len(self.ParametersArr)
        
        self.ParameterMinTxt = []
        self.ParameterValuesTxt = []
        self.ParameterMaxTxt = []
        self.Checkboxes = []
        self.ParameterNamesTxt = []
        self.ParameterSizer = []
        self.ParameterSpacers = []
        self.ParameterSpacersSizers = []
        
        for i in range(NumberOfParameters):
            self.ParameterMinTxt.append(0)
            self.ParameterValuesTxt.append(0)
            self.ParameterMaxTxt.append(0)
            self.Checkboxes.append(0)            
            self.ParameterNamesTxt.append(0)
            self.ParameterSizer.append(0)
            self.ParameterSpacers.append(0)
            self.ParameterSpacersSizers.append(0)
            
            self.ParameterSizer[i] = wx.BoxSizer(wx.HORIZONTAL)
            self.ParameterSizer[i].AddSpacer(10)
            
            self.ParameterMinTxt[i] = wx.TextCtrl(self.RightPanel, -1, '', size = (70, -1), style = wx.TE_CENTRE)
            self.ParameterSizer[i].Add(self.ParameterMinTxt[i])
            
            self.ParameterSizer[i].AddSpacer(20)
            
            self.ParameterValuesTxt[i] = wx.TextCtrl(self.RightPanel, -1, '', size = (100, -1), style = wx.TE_CENTRE)
            self.ParameterSizer[i].Add(self.ParameterValuesTxt[i])
            
            self.ParameterSizer[i].AddSpacer(20)
            
            self.ParameterMaxTxt[i] = wx.TextCtrl(self.RightPanel, -1, '', size = (70, -1), style = wx.TE_CENTRE)
            self.ParameterSizer[i].Add(self.ParameterMaxTxt[i])
            
            self.ParameterSizer[i].AddSpacer(20)
            
            self.Checkboxes[i] = wx.CheckBox(self.RightPanel, -1, '', size = (30, -1))
            self.ParameterSizer[i].Add(self.Checkboxes[i])
            
            self.ParameterSizer[i].AddSpacer(20)
            
            self.ParameterNamesTxt[i] = wx.StaticText(self.RightPanel, -1, '')
            self.ParameterSizer[i].Add(self.ParameterNamesTxt[i])
            
            self.RightBoxSizer.Add(self.ParameterSizer[i])
            
            self.ParameterSpacersSizers[i] = wx.BoxSizer(wx.VERTICAL)
            self.ParameterSpacersSizers[i].AddSpacer(3)
            
            self.ParameterSpacers[i] = wx.BoxSizer(wx.VERTICAL)
            self.ParameterSpacers[i].Add(self.ParameterSpacersSizers[i])
        
        self.RightPanel.Layout()
        self.Layout()
    
    def PrintParameters(self, filename):
        file = open(filename, 'rU')
        self.ParametersArr = file.readlines()
        file.close()
        
        self.ParameterMin = []
        self.ParameterValues = []
        self.ParameterMax = []
        self.ParameterOnOff = []
        self.ParameterNames = []
        
        for line in self.ParametersArr:
            LineArr = re.split(r'[\s]*', line)
            self.ParameterMin.append(LineArr[0])
            self.ParameterValues.append(LineArr[1])
            self.ParameterMax.append(LineArr[2])
            self.ParameterOnOff.append(LineArr[3])
            self.ParameterNames.append(LineArr[4])
        
        self.OnOff = []
        NumberOfParameters = len(self.ParameterValues)
        
        for i in range(NumberOfParameters):
            Font = wx.Font(9, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
            self.ParameterMinTxt[i].SetFont(Font)
            self.ParameterMinTxt[i].SetValue(self.ParameterMin[i])
            
            Font = wx.Font(11, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
            self.ParameterValuesTxt[i].SetFont(Font)
            self.ParameterValuesTxt[i].SetValue(self.ParameterValues[i])
            
            Font = wx.Font(9, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
            self.ParameterMaxTxt[i].SetFont(Font)
            self.ParameterMaxTxt[i].SetValue(self.ParameterMax[i])
            self.ParameterNamesTxt[i].SetLabel('%(#)02d: ' % {'#': i} + self.ParameterNames[i])
            
            if self.ParameterOnOff[i] == '1':
                self.OnOff.append(False)
            else:
                self.OnOff.append(True)
            
            self.Checkboxes[i].SetValue(self.OnOff[i])
    
    def RemoveParameters(self):
        for i in reversed(range(len(self.Checkboxes))):
            self.ParameterMinTxt[i].Destroy()
            self.ParameterValuesTxt[i].Destroy()
            self.ParameterMaxTxt[i].Destroy()
            self.Checkboxes[i].Destroy()
            self.ParameterNamesTxt[i].Destroy()
            self.RightBoxSizer.Remove(self.ParameterSizer[i])
            self.RightBoxSizer.Remove(self.ParameterSpacers[i])

## Function for copying results to appropriate directories
    def PrintResultsToCardfileDirectory(self):
        DataArr = []
        FitArr = []
        
        for root, dirs, files in os.walk('./'):
            for name in files:
                Filename = os.path.join(root, name)
                
                if os.path.splitext(Filename)[1] == '.dat':
                    ShortName = os.path.basename(Filename)
                    
                    if ShortName[0:3] == 'fit':
                        FitArr.append(os.path.join('.data/', ShortName))
                        
                    if ShortName[0:3] == 'dat':
                        DataArr.append(os.path.join('.data/', ShortName))
        
        for file in FitArr:
            try:
                shutil.copy(file, self.ResultsDirectory)
            except:
                pass
        
        for file in DataArr:
            try:
                shutil.copy(file, self.ResultsDirectory)
            except:
                pass
        
        try:
            shutil.copy(os.path.join('.data/', 'CovarianceMatrix.dat'), self.ResultsDirectory)
        except:
            pass
        
        try:
            shutil.copy('ParametersAfter.mcp', os.path.join(self.ResultsDirectory, 'LastFit.par'))
        except:
            pass
    
## Function for clearing directory of old data-files and parameter-files
    def ClearOldData(self):
        try:
            shutil.rmtree('./.data/')
        except:
            pass
    
    def ClearOldParameters(self):
        try:
            os.remove('ParametersBefore.mcp')
        except:
            pass
        
        try:
            os.remove('ParametersAfter.mcp')
        except:
            pass
 
## Define exit functions
    def CloseWindowFnc(self, event):        
        try:
            self.ClearOldData()
        except:
            pass
        
        try:
            self.ClearOldParameters()
        except:
            pass
        
        try:
            pylab.close()
        except:
            pass
        
        sys.exit(0)

class WaitingDialog(wx.Dialog):
    def __init__(self, Message, Title, Process):
        wx.Dialog.__init__(self, None, -1, Title, size=(300, 80))
        self.CenterOnScreen(wx.BOTH)
        
        self.RunningProcess = Process
        
        Infotext = wx.StaticText(self, -1, Message)
        CancelButton = wx.Button(self, wx.ID_CANCEL, "Abort")
        
        BoxSizer = wx.BoxSizer(wx.VERTICAL)
        BoxSizer.Add(Infotext, 1, wx.ALIGN_CENTER|wx.TOP, 10)
        BoxSizer.Add(CancelButton, 1, wx.ALIGN_CENTER|wx.BOTTOM, 10)
        
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.SetSizer(BoxSizer)
        
        self.Timer = wx.Timer(self)
        self.Timer.Start(1000)
        self.Bind(wx.EVT_TIMER, self.OnTimer, self.Timer)  

    def OnTimer(self, event):  
        if self.RunningProcess.poll() == 0:
            if self.IsModal():
                self.EndModal(wx.ID_OK)
            self.Timer.Stop()
            self.Destroy()

    def OnCancel(self):
        self.EndModal(wx.ID_CANCEL)
        self.Timer.Stop()
        self.Destroy()
    
    def OnClose(self, event):
        self.EndModal(wx.ID_CANCEL)
        self.Timer.Stop()
        self.Destroy()

## Read data files
def GetData(filename):
    file = open(filename, 'rU')
    data = file.readlines()
    file.close()
    
    DataStrArr = []

    for line in data:
        DataStrArr.append(line)

    FileExtensionStr = os.path.splitext(filename)[1]
    LineWithNumberOfDatapoint = DataStrArr[1]
    NumberOfDatapoints = int(re.split('=', LineWithNumberOfDatapoint)[1])

    xDataStr = []
    yDataStr = []
    wDataStr = []
    
    xDataFlt = []
    yDataFlt = []
    wDataFlt = []

    Contrast = float(re.split(r'=', DataStrArr[2])[1])
    
    for i in range(NumberOfDatapoints):
        xDataStr.append(0)
        yDataStr.append(0)
        wDataStr.append(0)

    for i in range(3, NumberOfDatapoints + 3):
        xDataStr[i - 3] = re.split(r'[\s]*', DataStrArr[i])[1]
        yDataStr[i - 3] = re.split(r'[\s]*', DataStrArr[i])[2]
        wDataStr[i - 3] = re.split(r'[\s]*', DataStrArr[i])[3]

    for item in xDataStr:
        xDataFlt.append(float(item))

    for item in yDataStr:
        yDataFlt.append(float(item))
        
    for item in wDataStr:
        wDataFlt.append(float(item))
    
    return xDataFlt, yDataFlt, wDataFlt, Contrast

def GetChisquare():
    file = open(os.path.join('./.data', 'ReturnMessage.mcp'), 'rU')
    Chisquare = file.readlines()[0]
    file.close()
    
    return Chisquare

def GetReturnMessage():
    file = open(os.path.join('./.data', 'ReturnMessage.mcp'), 'rU')
    ReturnMessage = file.readline()
    file.close()
    
    return ReturnMessage

## Functions for looking up plot colors
def FormatFnc(i):
    if i <= 100 and i > 80:
        return 'g.'
    elif i <= 80 and i > 60:
        return 'c.'
    elif i <= 60 and i > 40:
        return 'b.'
    elif i <= 40 and i > 20:
        return 'm.'
    elif i <= 20 and i >= 0:
        return 'y.'
    else:
        return 'r.'

## Boilerplate code linking program and widgets
if __name__ == '__main__':
    App = wx.App()
    Frame = MainCls(parent = None, id = -1)

    IconPath = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), 'WillItFit.png')
    Icon = wx.Icon()
    Icon.CopyFromBitmap(wx.Bitmap(IconPath, wx.BITMAP_TYPE_ANY))
    Frame.SetIcon(Icon)

    Frame.Show()
    App.MainLoop()
