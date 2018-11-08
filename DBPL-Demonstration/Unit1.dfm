object Form1: TForm1
  Left = 0
  Top = 0
  BorderStyle = bsDialog
  Caption = 'About'
  ClientHeight = 204
  ClientWidth = 368
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  Position = poMainFormCenter
  PixelsPerInch = 96
  TextHeight = 13
  object Panel1: TPanel
    Left = 0
    Top = 0
    Width = 368
    Height = 153
    Align = alTop
    BevelInner = bvLowered
    TabOrder = 0
    ExplicitWidth = 358
    object Label1: TLabel
      Left = 40
      Top = 28
      Width = 297
      Height = 41
      Caption = 
        'Biliniarity pairings Demo. Implemented library for several BN,BL' +
        'S,KSS and MNT curves.'
      WordWrap = True
    end
    object Label2: TLabel
      Left = 40
      Top = 70
      Width = 273
      Height = 104
      Caption = 
        'Implemented by FARAOUN Kamel Mohamed. EDDIS Laboratory. UDL-Univ' +
        'ersity- Algeria.'
      WordWrap = True
    end
    object Label3: TLabel
      Left = 152
      Top = 109
      Width = 24
      Height = 13
      Caption = '2018'
      WordWrap = True
    end
  end
  object Button1: TButton
    Left = 136
    Top = 167
    Width = 75
    Height = 25
    Caption = '&Ok'
    TabOrder = 1
    OnClick = Button1Click
  end
end
