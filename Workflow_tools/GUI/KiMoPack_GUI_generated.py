# -*- coding: utf-8 -*-

###########################################################################
## Python code generated with wxFormBuilder (version 3.10.1-0-g8feb16b3)
## http://www.wxformbuilder.org/
##
## PLEASE DO *NOT* EDIT THIS FILE!
###########################################################################

import wx
import wx.xrc
import wx.aui

###########################################################################
## Class frameMain
###########################################################################

class frameMain ( wx.Frame ):

	def __init__( self, parent ):
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = wx.EmptyString, pos = wx.DefaultPosition, size = wx.Size( -1,-1 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )

		self.SetSizeHints( wx.Size( -1,550 ), wx.DefaultSize )

		self.m_toolBar1 = self.CreateToolBar( wx.TB_HORIZONTAL, wx.ID_ANY )
		self.working_directory_label = wx.StaticText( self.m_toolBar1, wx.ID_ANY, u"current working directory", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.working_directory_label.Wrap( -1 )

		self.m_toolBar1.AddControl( self.working_directory_label )
		self.working_directory = wx.DirPickerCtrl( self.m_toolBar1, wx.ID_ANY, wx.EmptyString, u"Select a folder", wx.DefaultPosition, wx.Size( 500,-1 ), wx.DIRP_DEFAULT_STYLE )
		self.working_directory.SetMinSize( wx.Size( 300,-1 ) )

		self.m_toolBar1.AddControl( self.working_directory )
		self.m_toolBar1.AddSeparator()

		self.m_toolBar1.AddSeparator()

		self.load_setting_label = wx.StaticText( self.m_toolBar1, wx.ID_ANY, u"Load setting", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.load_setting_label.Wrap( -1 )

		self.m_toolBar1.AddControl( self.load_setting_label )
		self.Load_settings = self.m_toolBar1.AddTool( wx.ID_ANY, u"tool", wx.Bitmap( u"embedded_files/document-open.png", wx.BITMAP_TYPE_ANY ), wx.NullBitmap, wx.ITEM_NORMAL, wx.EmptyString, wx.EmptyString, None )

		self.m_toolBar1.AddSeparator()

		self.m_toolBar1.AddSeparator()

		self.save_setting_label = wx.StaticText( self.m_toolBar1, wx.ID_ANY, u"save settings", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.save_setting_label.Wrap( -1 )

		self.m_toolBar1.AddControl( self.save_setting_label )
		self.save_settings = self.m_toolBar1.AddTool( wx.ID_ANY, u"tool", wx.Bitmap( u"embedded_files/document-save.png", wx.BITMAP_TYPE_ANY ), wx.NullBitmap, wx.ITEM_NORMAL, u"save current settings", wx.EmptyString, None )

		self.m_toolBar1.AddSeparator()

		self.m_toolBar1.AddSeparator()

		self.send_setting_label = wx.StaticText( self.m_toolBar1, wx.ID_ANY, u"send settings", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.send_setting_label.Wrap( -1 )

		self.m_toolBar1.AddControl( self.send_setting_label )
		self.send_email_to_creator = self.m_toolBar1.AddTool( wx.ID_ANY, u"tool", wx.Bitmap( u"embedded_files/emblem-mail.png", wx.BITMAP_TYPE_ANY ), wx.NullBitmap, wx.ITEM_NORMAL, u"dump settings in email and ask for help", wx.EmptyString, None )

		self.m_toolBar1.Realize()

		frameMain_sizer = wx.BoxSizer( wx.HORIZONTAL )

		wSizer1 = wx.WrapSizer( wx.HORIZONTAL, wx.WRAPSIZER_DEFAULT_FLAGS )

		self.frameMainNotebook = wx.aui.AuiNotebook( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.aui.AUI_NB_DEFAULT_STYLE )
		self.data_loading_panel = wx.Panel( self.frameMainNotebook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		gbSizer2 = wx.GridBagSizer( 10, 10 )
		gbSizer2.SetFlexibleDirection( wx.BOTH )
		gbSizer2.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )

		self.open_single_file_gui_button = wx.Button( self.data_loading_panel, wx.ID_ANY, u"Load data from single file with GUI", wx.Point( -1,-1 ), wx.DefaultSize, 0 )
		self.open_single_file_gui_button.SetBackgroundColour( wx.SystemSettings.GetColour( wx.SYS_COLOUR_HIGHLIGHT ) )
		self.open_single_file_gui_button.SetMinSize( wx.Size( -1,40 ) )

		gbSizer2.Add( self.open_single_file_gui_button, wx.GBPosition( 0, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )

		self.Load_recent = wx.Button( self.data_loading_panel, wx.ID_ANY, u"Reload last used file", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.Load_recent.SetBackgroundColour( wx.SystemSettings.GetColour( wx.SYS_COLOUR_HIGHLIGHT ) )
		self.Load_recent.SetMinSize( wx.Size( -1,40 ) )

		gbSizer2.Add( self.Load_recent, wx.GBPosition( 1, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )

		self.load_filename1 = wx.Button( self.data_loading_panel, wx.ID_ANY, u"Load file with name", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.load_filename1.SetBackgroundColour( wx.SystemSettings.GetColour( wx.SYS_COLOUR_HIGHLIGHT ) )
		self.load_filename1.SetMinSize( wx.Size( -1,40 ) )

		gbSizer2.Add( self.load_filename1, wx.GBPosition( 2, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )

		self.load_filename = wx.TextCtrl( self.data_loading_panel, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.Size( 200,-1 ), 0 )
		gbSizer2.Add( self.load_filename, wx.GBPosition( 3, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )

		bSizer2 = wx.BoxSizer( wx.VERTICAL )

		type_SIAChoices = [ u"SIA", u"hdf5", u"custom" ]
		self.type_SIA = wx.RadioBox( self.data_loading_panel, wx.ID_ANY, u"Filetype", wx.DefaultPosition, wx.DefaultSize, type_SIAChoices, 1, wx.RA_SPECIFY_COLS )
		self.type_SIA.SetSelection( 0 )
		bSizer2.Add( self.type_SIA, 0, wx.ALL, 5 )

		self.custom_filetype = wx.TextCtrl( self.data_loading_panel, wx.ID_ANY, u"SIA", wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer2.Add( self.custom_filetype, 0, wx.ALL, 5 )


		gbSizer2.Add( bSizer2, wx.GBPosition( 4, 0 ), wx.GBSpan( 3, 1 ), wx.EXPAND, 5 )

		self.path_label = wx.StaticText( self.data_loading_panel, wx.ID_ANY, u"From path", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.path_label.Wrap( -1 )

		gbSizer2.Add( self.path_label, wx.GBPosition( 7, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )

		self.load_path = wx.TextCtrl( self.data_loading_panel, wx.ID_ANY, u"None", wx.DefaultPosition, wx.Size( 200,-1 ), 0 )
		gbSizer2.Add( self.load_path, wx.GBPosition( 8, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )

		loading_option_sizer = wx.GridSizer( 11, 2, 0, 0 )

		self.Data_separator_label = wx.StaticText( self.data_loading_panel, wx.ID_ANY, u"Data separator", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.Data_separator_label.Wrap( -1 )

		loading_option_sizer.Add( self.Data_separator_label, 0, wx.ALL, 5 )

		self.data_separator = wx.TextCtrl( self.data_loading_panel, wx.ID_ANY, u"\\t", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.data_separator.SetMaxSize( wx.Size( 30,-1 ) )

		loading_option_sizer.Add( self.data_separator, 0, wx.ALL, 5 )

		self.data_decimal_separator_label = wx.StaticText( self.data_loading_panel, wx.ID_ANY, u"decimal separator", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.data_decimal_separator_label.Wrap( -1 )

		loading_option_sizer.Add( self.data_decimal_separator_label, 0, wx.ALL, 5 )

		self.data_decimal_separator = wx.TextCtrl( self.data_loading_panel, wx.ID_ANY, u",", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.data_decimal_separator.SetMaxSize( wx.Size( 30,-1 ) )

		loading_option_sizer.Add( self.data_decimal_separator, 0, wx.ALL, 5 )

		self.index_is_energy_label = wx.StaticText( self.data_loading_panel, wx.ID_ANY, u"Index is Energy", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.index_is_energy_label.Wrap( -1 )

		loading_option_sizer.Add( self.index_is_energy_label, 0, wx.ALL, 5 )

		self.Index_is_energy_check = wx.CheckBox( self.data_loading_panel, wx.ID_ANY, u"convert to nm", wx.DefaultPosition, wx.DefaultSize, 0 )
		loading_option_sizer.Add( self.Index_is_energy_check, 0, wx.ALL, 5 )

		self.transpose_label = wx.StaticText( self.data_loading_panel, wx.ID_ANY, u"Transpose", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.transpose_label.Wrap( -1 )

		loading_option_sizer.Add( self.transpose_label, 0, wx.ALL, 5 )

		self.transpose_checkbox = wx.CheckBox( self.data_loading_panel, wx.ID_ANY, u"swap axis", wx.DefaultPosition, wx.DefaultSize, 0 )
		loading_option_sizer.Add( self.transpose_checkbox, 0, wx.ALL, 5 )

		self.sort_indexes_label = wx.StaticText( self.data_loading_panel, wx.ID_ANY, u"Resort the indexes", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.sort_indexes_label.Wrap( -1 )

		loading_option_sizer.Add( self.sort_indexes_label, 0, wx.ALL, 5 )

		self.resort_indexes_check = wx.CheckBox( self.data_loading_panel, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, 0 )
		loading_option_sizer.Add( self.resort_indexes_check, 0, wx.ALL, 5 )

		self.divide_times_label = wx.StaticText( self.data_loading_panel, wx.ID_ANY, u"divide times by", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.divide_times_label.Wrap( -1 )

		loading_option_sizer.Add( self.divide_times_label, 0, wx.ALL, 5 )

		self.shift_times_by = wx.TextCtrl( self.data_loading_panel, wx.ID_ANY, u"None", wx.DefaultPosition, wx.Size( 50,-1 ), 0 )
		self.shift_times_by.SetMinSize( wx.Size( 110,-1 ) )

		loading_option_sizer.Add( self.shift_times_by, 0, wx.ALL, 5 )

		self.shift_times_label = wx.StaticText( self.data_loading_panel, wx.ID_ANY, u"shift times by", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.shift_times_label.Wrap( -1 )

		loading_option_sizer.Add( self.shift_times_label, 0, wx.ALL, 5 )

		self.shift_times_value = wx.TextCtrl( self.data_loading_panel, wx.ID_ANY, u"None", wx.DefaultPosition, wx.Size( 60,-1 ), 0 )
		self.shift_times_value.SetMinSize( wx.Size( 110,-1 ) )

		loading_option_sizer.Add( self.shift_times_value, 0, wx.ALL, 5 )

		self.data_type_label = wx.StaticText( self.data_loading_panel, wx.ID_ANY, u"intensity axis [Delta OD]", wx.DefaultPosition, wx.Size( -1,-1 ), 0 )
		self.data_type_label.Wrap( -1 )

		loading_option_sizer.Add( self.data_type_label, 0, wx.ALL, 5 )

		self.intensity_axis_values = wx.TextCtrl( self.data_loading_panel, wx.ID_ANY, u"None", wx.DefaultPosition, wx.Size( -1,-1 ), 0 )
		self.intensity_axis_values.SetMinSize( wx.Size( 110,-1 ) )

		loading_option_sizer.Add( self.intensity_axis_values, 0, wx.ALL, 5 )

		self.time_units_label = wx.StaticText( self.data_loading_panel, wx.ID_ANY, u"time units [ps] default", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.time_units_label.Wrap( -1 )

		loading_option_sizer.Add( self.time_units_label, 0, wx.ALL, 5 )

		self.time_units = wx.TextCtrl( self.data_loading_panel, wx.ID_ANY, u"None", wx.DefaultPosition, wx.DefaultSize, 0 )
		loading_option_sizer.Add( self.time_units, 0, wx.ALL, 5 )

		self.external_time_label = wx.StaticText( self.data_loading_panel, wx.ID_ANY, u"external time file", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.external_time_label.Wrap( -1 )

		loading_option_sizer.Add( self.external_time_label, 0, wx.ALL, 5 )

		self.external_time_file = wx.TextCtrl( self.data_loading_panel, wx.ID_ANY, u"None", wx.DefaultPosition, wx.DefaultSize, 0 )
		loading_option_sizer.Add( self.external_time_file, 0, wx.ALL, 5 )

		self.external_wave_file_label = wx.StaticText( self.data_loading_panel, wx.ID_ANY, u"external wavelength file", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.external_wave_file_label.Wrap( -1 )

		loading_option_sizer.Add( self.external_wave_file_label, 0, wx.ALL, 5 )

		self.external_wavelength_file = wx.TextCtrl( self.data_loading_panel, wx.ID_ANY, u"None", wx.DefaultPosition, wx.DefaultSize, 0 )
		loading_option_sizer.Add( self.external_wavelength_file, 0, wx.ALL, 5 )


		gbSizer2.Add( loading_option_sizer, wx.GBPosition( 0, 1 ), wx.GBSpan( 8, 2 ), wx.EXPAND, 5 )


		gbSizer2.AddGrowableCol( 2 )
		gbSizer2.AddGrowableRow( 8 )

		self.data_loading_panel.SetSizer( gbSizer2 )
		self.data_loading_panel.Layout()
		gbSizer2.Fit( self.data_loading_panel )
		self.frameMainNotebook.AddPage( self.data_loading_panel, u"Data Loading", True, wx.NullBitmap )
		self.data_cleaning_panel = wx.Panel( self.frameMainNotebook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.frameMainNotebook.AddPage( self.data_cleaning_panel, u"Data Cleaning", False, wx.NullBitmap )
		self.RAW_plotting_panel = wx.Panel( self.frameMainNotebook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.frameMainNotebook.AddPage( self.RAW_plotting_panel, u"RAW plotting", False, wx.NullBitmap )
		self.data_fitting_panel = wx.Panel( self.frameMainNotebook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.frameMainNotebook.AddPage( self.data_fitting_panel, u"Fitting", False, wx.NullBitmap )
		self.Single_scan_panel = wx.Panel( self.frameMainNotebook, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.frameMainNotebook.AddPage( self.Single_scan_panel, u"Single Scan handling", False, wx.NullBitmap )

		wSizer1.Add( self.frameMainNotebook, 1, wx.EXPAND |wx.ALL, 5 )

		bSizer3 = wx.BoxSizer( wx.HORIZONTAL )

		self.Logo = wx.StaticBitmap( self, wx.ID_ANY, wx.Bitmap( u"embedded_files/KiMoPack_logo.png", wx.BITMAP_TYPE_ANY ), wx.DefaultPosition, wx.DefaultSize, 0 )
		bSizer3.Add( self.Logo, 0, wx.ALL, 5 )

		self.copyright = wx.StaticText( self, wx.ID_ANY, u"copyright by Jens Uhlig 2022 \n Jens.uhlig@chemphys.lu.se \n www.chemphys.lu.se/research/~ \n /projects/kimopack/", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.copyright.Wrap( -1 )

		bSizer3.Add( self.copyright, 0, wx.ALL, 5 )

		self.QR = wx.StaticBitmap( self, wx.ID_ANY, wx.Bitmap( u"embedded_files/qr-code50.png", wx.BITMAP_TYPE_ANY ), wx.DefaultPosition, wx.Size( -1,-1 ), 0 )
		self.QR.Hide()
		self.QR.SetMaxSize( wx.Size( 40,40 ) )

		bSizer3.Add( self.QR, 0, wx.ALL, 5 )


		wSizer1.Add( bSizer3, 1, wx.EXPAND, 5 )

		self.m_filePicker1 = wx.FilePickerCtrl( self, wx.ID_ANY, wx.EmptyString, u"Select a file", u"*.*", wx.DefaultPosition, wx.DefaultSize, wx.FLP_DEFAULT_STYLE )
		wSizer1.Add( self.m_filePicker1, 0, wx.ALL, 5 )


		frameMain_sizer.Add( wSizer1, 1, wx.EXPAND, 5 )

		self.log_screen = wx.TextCtrl( self, wx.ID_ANY, wx.EmptyString, wx.DefaultPosition, wx.DefaultSize, wx.HSCROLL|wx.TE_MULTILINE|wx.TE_READONLY|wx.HSCROLL|wx.VSCROLL )
		self.log_screen.SetMinSize( wx.Size( 500,550 ) )

		frameMain_sizer.Add( self.log_screen, 0, wx.ALL, 5 )


		self.SetSizer( frameMain_sizer )
		self.Layout()
		frameMain_sizer.Fit( self )

		self.Centre( wx.BOTH )

		# Connect Events
		self.working_directory.Bind( wx.EVT_DIRPICKER_CHANGED, self.working_directoryOnDirChanged )
		self.Bind( wx.EVT_TOOL, self.Load_settingsOnToolClicked, id = self.Load_settings.GetId() )
		self.Bind( wx.EVT_TOOL, self.save_settingsOnToolClicked, id = self.save_settings.GetId() )
		self.Bind( wx.EVT_TOOL, self.send_email_to_creatorOnToolClicked, id = self.send_email_to_creator.GetId() )
		self.open_single_file_gui_button.Bind( wx.EVT_BUTTON, self.open_single_file_gui_buttonOnButtonClick )
		self.Load_recent.Bind( wx.EVT_BUTTON, self.Load_recentOnButtonClick )
		self.load_filename1.Bind( wx.EVT_BUTTON, self.load_filenameOnButtonClick )
		self.type_SIA.Bind( wx.EVT_RADIOBOX, self.type_SIAOnRadioBox )

	def __del__( self ):
		pass


	# Virtual event handlers, override them in your derived class
	def working_directoryOnDirChanged( self, event ):
		event.Skip()

	def Load_settingsOnToolClicked( self, event ):
		event.Skip()

	def save_settingsOnToolClicked( self, event ):
		event.Skip()

	def send_email_to_creatorOnToolClicked( self, event ):
		event.Skip()

	def open_single_file_gui_buttonOnButtonClick( self, event ):
		event.Skip()

	def Load_recentOnButtonClick( self, event ):
		event.Skip()

	def load_filenameOnButtonClick( self, event ):
		event.Skip()

	def type_SIAOnRadioBox( self, event ):
		event.Skip()


