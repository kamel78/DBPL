object Form4: TForm4
  Left = 0
  Top = 0
  BorderIcons = []
  BorderStyle = bsToolWindow
  Caption = 'Save As ...'
  ClientHeight = 150
  ClientWidth = 282
  Color = clWindow
  Ctl3D = False
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  Position = poMainFormCenter
  PixelsPerInch = 96
  TextHeight = 13
  object Button1: TButton
    Left = 32
    Top = 115
    Width = 75
    Height = 25
    Caption = '&Ok'
    Default = True
    TabOrder = 0
    OnClick = Button1Click
  end
  object Button2: TButton
    Left = 160
    Top = 115
    Width = 75
    Height = 25
    Cancel = True
    Caption = '&Cancel'
    TabOrder = 1
    OnClick = Button2Click
  end
  object GroupBox1: TGroupBox
    Left = 0
    Top = 0
    Width = 282
    Height = 105
    Align = alTop
    Caption = 'Output Format  '
    Ctl3D = True
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -11
    Font.Name = 'Tahoma'
    Font.Style = [fsBold]
    ParentCtl3D = False
    ParentFont = False
    TabOrder = 2
    object RadioButton2: TRadioButton
      Left = 16
      Top = 25
      Width = 265
      Height = 17
      Caption = 'Binary format (Shorter)'
      Checked = True
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -11
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
      TabOrder = 0
      TabStop = True
    end
    object RadioButton3: TRadioButton
      Left = 16
      Top = 50
      Width = 265
      Height = 17
      Caption = 'Base-64 Format (Text-compatible)'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -11
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
      TabOrder = 1
    end
    object CheckBox1: TCheckBox
      Left = 10
      Top = 77
      Width = 168
      Height = 17
      Caption = 'Include Secrete Key'
      Checked = True
      State = cbChecked
      TabOrder = 2
    end
  end
end
