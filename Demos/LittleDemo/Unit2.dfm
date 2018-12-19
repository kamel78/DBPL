object Form2: TForm2
  Left = 0
  Top = 0
  Caption = 'Form2'
  ClientHeight = 247
  ClientWidth = 1001
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Label1: TLabel
    Left = 16
    Top = 24
    Width = 106
    Height = 13
    Caption = 'First Identity (Point P)'
  end
  object Label2: TLabel
    Left = 15
    Top = 51
    Width = 122
    Height = 13
    Caption = 'Second Identity (Point Q)'
  end
  object Label3: TLabel
    Left = 15
    Top = 85
    Width = 95
    Height = 13
    Caption = 'Random number "a"'
  end
  object Label4: TLabel
    Left = 15
    Top = 109
    Width = 95
    Height = 13
    Caption = 'Random number "b"'
  end
  object Label5: TLabel
    Left = 15
    Top = 156
    Width = 101
    Height = 13
    Caption = 'Pairing e(P,Q)^(a*b)'
  end
  object Label6: TLabel
    Left = 16
    Top = 181
    Width = 91
    Height = 13
    Caption = 'Pairing e(a*P,b*Q)'
  end
  object Label7: TLabel
    Left = 15
    Top = 131
    Width = 72
    Height = 13
    Caption = 'Pairings e(P,Q)'
  end
  object Label8: TLabel
    Left = 120
    Top = 208
    Width = 28
    Height = 13
    Caption = 'Time'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -11
    Font.Name = 'Tahoma'
    Font.Style = [fsBold]
    ParentFont = False
  end
  object Button1: TButton
    Left = 298
    Top = 24
    Width = 121
    Height = 25
    Caption = 'Hash then Paire'
    TabOrder = 0
    OnClick = Button1Click
  end
  object Edit1: TEdit
    Left = 116
    Top = 154
    Width = 877
    Height = 21
    TabOrder = 1
  end
  object Edit2: TEdit
    Left = 115
    Top = 178
    Width = 878
    Height = 21
    TabOrder = 2
  end
  object Edit3: TEdit
    Left = 154
    Top = 21
    Width = 121
    Height = 21
    TabOrder = 3
    Text = 'Identity 1'
  end
  object Edit4: TEdit
    Left = 154
    Top = 48
    Width = 121
    Height = 21
    TabOrder = 4
    Text = 'Identity 2'
  end
  object Edit5: TEdit
    Left = 116
    Top = 82
    Width = 665
    Height = 21
    TabOrder = 5
  end
  object Edit6: TEdit
    Left = 116
    Top = 105
    Width = 665
    Height = 21
    TabOrder = 6
  end
  object Edit7: TEdit
    Left = 116
    Top = 129
    Width = 877
    Height = 21
    TabOrder = 7
  end
  object BNCurvePairing1: TBNCurvePairing
    CoordinatesSystem = csProjective
    Fp12PoweringMode = pmKarbina
    PairingAlgorithm = bpOptAte
    LoopMode = lpmAuto
    Parametres = scBCMNPZ
    Left = 464
    Top = 32
  end
end
