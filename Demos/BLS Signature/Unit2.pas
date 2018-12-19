unit Unit2;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.ExtCtrls, Vcl.ComCtrls,
  Vcl.StdCtrls, Vcl.Buttons, BNCurves, BLS24Curves, GeneralTypes, BLS12Curves,
  KSS36Curves, BLS27Curves, KSS18Curves, SuperSingularCurves, MNTCurves,
  KSS16Curves,ECCFp,ECCFp2,ECCFp3,ECCFp4,ECCFp6,ECCFp9, LargeIntegers, System.IOUtils, Hashfunctions;

type
  TForm2 = class(TForm)
    Panel1: TPanel;
    Bevel1: TBevel;
    Label31: TLabel;
    Label32: TLabel;
    Label39: TLabel;
    Label40: TLabel;
    Label41: TLabel;
    SpeedButton5: TSpeedButton;
    SpeedButton6: TSpeedButton;
    Label44: TLabel;
    Label45: TLabel;
    Label46: TLabel;
    Button1: TButton;
    Edit22: TEdit;
    Edit23: TEdit;
    Edit30: TEdit;
    Edit31: TEdit;
    ComboBox10: TComboBox;
    ComboBox11: TComboBox;
    Label42: TLabel;
    Bevel4: TBevel;
    Label1: TLabel;
    ComboBox1: TComboBox;
    PageControl1: TPageControl;
    TabSheet1: TTabSheet;
    TabSheet2: TTabSheet;
    Label33: TLabel;
    Label34: TLabel;
    Label35: TLabel;
    Label36: TLabel;
    Label43: TLabel;
    SpeedButton1: TSpeedButton;
    SpeedButton2: TSpeedButton;
    Image1: TImage;
    SpeedButton3: TSpeedButton;
    SpeedButton4: TSpeedButton;
    Button2: TButton;
    Button3: TButton;
    Edit24: TEdit;
    Edit25: TEdit;
    Edit26: TEdit;
    Edit27: TEdit;
    Edit32: TEdit;
    Bevel3: TBevel;
    Bevel2: TBevel;
    BLS12CurvePairing1: TBLS12CurvePairing;
    BLS24CurvePairing1: TBLS24CurvePairing;
    BNCurvePairing1: TBNCurvePairing;
    KSS16CurvePairing1: TKSS16CurvePairing;
    MNTCurve1: TMNTCurve;
    SuperSingularCurvePairing1: TSuperSingularCurvePairing;
    KSS18CurvePairing1: TKSS18CurvePairing;
    BLS27Curve1: TBLS27Curve;
    KSS36Curve1: TKSS36Curve;
    OpenDialog1: TOpenDialog;
    SaveDialog1: TSaveDialog;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Memo1: TMemo;
    Edit1: TEdit;
    Button4: TButton;
    Button5: TButton;
    Label6: TLabel;
    Label7: TLabel;
    Label8: TLabel;
    Label9: TLabel;
    Label10: TLabel;
    Label11: TLabel;
    procedure ComboBox10Change(Sender: TObject);
    procedure FormShow(Sender: TObject);
    procedure Button1Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
    procedure SpeedButton1Click(Sender: TObject);
    procedure SpeedButton2Click(Sender: TObject);
    procedure Button3Click(Sender: TObject);
    procedure SpeedButton3Click(Sender: TObject);
    procedure SpeedButton4Click(Sender: TObject);
    procedure ComboBox1Change(Sender: TObject);
    procedure ComboBox11Change(Sender: TObject);
    procedure SpeedButton6Click(Sender: TObject);
    procedure SpeedButton5Click(Sender: TObject);
    procedure Button4Click(Sender: TObject);
    procedure Button5Click(Sender: TObject);
  private
    { Déclarations privées }
  public
    { Déclarations publiques }
  end;

var
  Form2: TForm2;

  G2_1,PubKeyInEFp:FpPoint;
  G2_2,PubKeyInEFp2:Fp2Point;
  G2_3,PubKeyInEFp3:Fp3Point;
  G2_4,PubKeyInEFp4:Fp4Point;
  G2_6,PubKeyInEFp6:Fp6Point;
  G2_9,PubKeyInEFp9:Fp9Point;
  SecretKey:LargeInt;
  e1_12,e2_12:GTBLS12;
  e1_16,e2_16:GTKSS16;
  e1_18,e2_18:GTKSS18;
  e1_27,e2_27:GTBLS27;
  e1_36,e2_36:GTKSS36;
  e1_24,e2_24:GtBls24;
  e1_6,e2_6:GTMNT;
  e1_2,e2_2:GtSS;

implementation

{$R *.dfm}

procedure TForm2.Button1Click(Sender: TObject);
var G1Size,Pblength:integer;
begin
case ComboBox10.ItemIndex of
0:begin // Super Singular Curves
  PubKeyInEFp.SetCurveParams(SuperSingularCurvePairing1.CurveParams);
  PubKeyInEFp:=SuperSingularCurvePairing1.GetDefautGGenerator;
  Edit30.Text:=PubKeyInEFp.ToHexString;
  SecretKey.SetToRandom(SuperSingularCurvePairing1.CurveParams.R);
  PubKeyInEFp:=SecretKey*PubKeyInEFp;
  Edit23.Text:=PubKeyInEFp.ToHexString;
  Edit31.Text:=ArrayToHex(PubKeyInEFp.CompressToArray);
  G1Size:=SuperSingularCurvePairing1.CurveParams.P.BitLength;
  Pblength:=PubKeyInEFp.X.BitLength;
  end;
1:begin   // MNT Curves
  PubKeyInEFp3.SetCurveParams(MNTCurve1.CurveParams);
  PubKeyInEFp3.SetToDefaultGenerator;
  Edit30.Text:=PubKeyInEFp3.ToHexString;
  SecretKey.SetToRandom(MNTCurve1.CurveParams.R);
  PubKeyInEFp3:=SecretKey*PubKeyInEFp3;
  Edit23.Text:=PubKeyInEFp3.ToHexString;
  Edit31.Text:=ArrayToHex(PubKeyInEFp3.CompressToArray);
  G1Size:=MNTCurve1.CurveParams.P.BitLength;
  Pblength:=PubKeyInEFp3.X.a.BitLength*3;
  end;
2:begin   // BN Curves
  PubKeyInEFp2.SetCurveParams(BNCurvePairing1.CurveParams);
  PubKeyInEFp2.SetToDefaultGenerator;
  Edit30.Text:=PubKeyInEFp2.ToHexString;
  SecretKey.SetToRandom(BNCurvePairing1.CurveParams.R);
  PubKeyInEFp2:=SecretKey*PubKeyInEFp2;
  Edit23.Text:=PubKeyInEFp2.ToHexString;
  Edit31.Text:=ArrayToHex(PubKeyInEFp2.CompressToArray);
  G1Size:=BNCurvePairing1.CurveParams.P.BitLength;
  Pblength:=PubKeyInEFp2.X.a.BitLength*2;
  end;
3:begin  // BLS12 Curves
  PubKeyInEFp2.SetCurveParams(BLS12CurvePairing1.CurveParams);
  PubKeyInEFp2.SetToDefaultGenerator;
  Edit30.Text:=PubKeyInEFp2.ToHexString;
  SecretKey.SetToRandom(BLS12CurvePairing1.CurveParams.R);
  PubKeyInEFp2:=SecretKey*PubKeyInEFp2;
  Edit23.Text:=PubKeyInEFp2.ToHexString;
  Edit31.Text:=ArrayToHex(PubKeyInEFp2.CompressToArray);
  G1Size:=BLS12CurvePairing1.CurveParams.P.BitLength;
  Pblength:=PubKeyInEFp2.X.a.BitLength*2;
  end;
4:begin   // BLS24 Curves
  PubKeyInEFp4.SetCurveParams(BLS24CurvePairing1.CurveParams);
  PubKeyInEFp4.SetToDefaultGenerator;
  Edit30.Text:=PubKeyInEFp4.ToHexString;
  SecretKey.SetToRandom(BLS24CurvePairing1.CurveParams.R);
  PubKeyInEFp4:=SecretKey*PubKeyInEFp4;
  Edit23.Text:=PubKeyInEFp4.ToHexString;
  Edit31.Text:=ArrayToHex(PubKeyInEFp4.CompressToArray);
  G1Size:=BLS24CurvePairing1.CurveParams.P.BitLength;
  Pblength:=PubKeyInEFp4.X.a.a.BitLength*4;
  end;
5:begin  // BLS27 Curves
  PubKeyInEFp9.SetCurveParams(BLS27Curve1.CurveParams);
  PubKeyInEFp9.SetToDefaultGenerator;
  Edit30.Text:=PubKeyInEFp9.ToHexString;
  SecretKey.SetToRandom(BLS27Curve1.CurveParams.R);
  PubKeyInEFp9:=SecretKey*PubKeyInEFp9;
  Edit23.Text:=PubKeyInEFp9.ToHexString;
  Edit31.Text:=ArrayToHex(PubKeyInEFp9.CompressToArray);
  G1Size:=BLS27Curve1.CurveParams.P.BitLength;
  Pblength:=PubKeyInEFp9.X.a.a.BitLength*9;
  end;
6:begin  // KSS16 Curves
  PubKeyInEFp4.SetCurveParams(KSS16CurvePairing1.CurveParams);
  PubKeyInEFp4.SetToDefaultGenerator;
  Edit30.Text:=PubKeyInEFp4.ToHexString;
  SecretKey.SetToRandom(KSS16CurvePairing1.CurveParams.R);
  PubKeyInEFp4:=SecretKey*PubKeyInEFp4;
  Edit23.Text:=PubKeyInEFp4.ToHexString;
  Edit31.Text:=ArrayToHex(PubKeyInEFp4.CompressToArray);
  G1Size:=KSS16CurvePairing1.CurveParams.P.BitLength;
  Pblength:=PubKeyInEFp4.X.a.a.BitLength*4;
  end;
7:begin // KSS18 Curves
  PubKeyInEFp3.SetCurveParams(KSS18CurvePairing1.CurveParams);
  PubKeyInEFp3.SetToDefaultGenerator;
  Edit30.Text:=PubKeyInEFp3.ToHexString;
  SecretKey.SetToRandom(KSS18CurvePairing1.CurveParams.R);
  PubKeyInEFp3:=SecretKey*PubKeyInEFp3;
  Edit23.Text:=PubKeyInEFp3.ToHexString;
  Edit31.Text:=ArrayToHex(PubKeyInEFp3.CompressToArray);
  G1Size:=KSS18CurvePairing1.CurveParams.P.BitLength;
  Pblength:=PubKeyInEFp3.X.a.BitLength*3;
  end;
8:begin  // KSS36 Curves
  PubKeyInEFp6.SetCurveParams(KSS36Curve1.CurveParams);
  PubKeyInEFp6.SetToDefaultGenerator;
  Edit30.Text:=PubKeyInEFp6.ToHexString;
  SecretKey.SetToRandom(KSS36Curve1.CurveParams.R);
  PubKeyInEFp6:=SecretKey*PubKeyInEFp6;
  Edit23.Text:=PubKeyInEFp6.ToHexString;
  Edit31.Text:=ArrayToHex(PubKeyInEFp6.CompressToArray);
  G1Size:=KSS36Curve1.CurveParams.P.BitLength;
  Pblength:=PubKeyInEFp6.X.a.a.BitLength*6;
  end;
end;
label31.Caption:='Secrete Key (an element from FP) : on '+inttostr(SecretKey.BitLength)+' bit';
Edit22.Text:=SecretKey.ToHexString;
Button2.Enabled:=true;
Button3.Enabled:=true;
Button4.Enabled:=true;
Button5.Enabled:=true;
SpeedButton5.Enabled:=true;
SpeedButton6.Enabled:=true;
label7.Caption:=inttostr(G1Size)+' bit';
label10.Caption:=inttostr(Pblength)+' bit';
label9.Caption:=inttostr(G1Size div 2)+' bit';
label8.Caption:=inttostr(G1Size)+' bit';
end;

procedure TForm2.Button2Click(Sender: TObject);
var FileHash,FileSignature:TBytes;
    SigAsPoint:FpPoint;
begin
if FileExists(Edit24.Text) then begin
                                if Tpath.HasValidPathChars(edit27.Text,true) then
                                                begin
                                                FileHash:=SHA256FileHash(Edit24.Text);
                                                case ComboBox10.ItemIndex of
                                                0:SigAsPoint.SetCurveParams(SuperSingularCurvePairing1.CurveParams);
                                                1:SigAsPoint.SetCurveParams(MNTCurve1.CurveParams);
                                                2:SigAsPoint.SetCurveParams(BNCurvePairing1.CurveParams);
                                                3:SigAsPoint.SetCurveParams(BLS12CurvePairing1.CurveParams);
                                                4:SigAsPoint.SetCurveParams(BLS24CurvePairing1.CurveParams);
                                                5:SigAsPoint.SetCurveParams(BLS27Curve1.CurveParams);
                                                6:SigAsPoint.SetCurveParams(KSS16CurvePairing1.CurveParams);
                                                7:SigAsPoint.SetCurveParams(KSS18CurvePairing1.CurveParams);
                                                8:SigAsPoint.SetCurveParams(KSS36Curve1.CurveParams);
                                                end;
                                                SigAsPoint.SetAsTorsionFromHash(FileHash);
                                                SigAsPoint:=SecretKey*SigAsPoint;
                                                FileSignature:=SigAsPoint.CompressToArray;
                                                TFile.WriteAllBytes(Edit27.Text,FileSignature);
                                                Edit32.Text:=ArrayToHex(FileSignature);
                                                Label46.Caption:='('+inttostr((length(FileSignature)-1)*8)+' bit)';
                                                MessageDlg('Signature Générée avec succés !',mtConfirmation,[mbok],0);
                                                end
                                else MessageDlg('Fichier destination de la signature invalide ',mtwarning,[mbok],0);
                                end
else MessageDlg('Fichier source à signé introuvable ',mtwarning,[mbok],0);
end;

procedure TForm2.Button3Click(Sender: TObject);
var filename:string;
    FileSignature,Filehash:TBytes;
    SigAsPoint:FpPoint;
    G2:Fp2Point;
    hashAspoint:FpPoint;
    Equals:boolean;
begin
if FileExists(Edit25.Text) then begin
                                filename:= Edit25.Text;
                                if FileExists(Edit26.Text) then begin
                                                case ComboBox10.ItemIndex of
                                                0:begin
                                                    SigAsPoint.SetCurveParams(SuperSingularCurvePairing1.CurveParams);
                                                    hashAspoint.SetCurveParams(SuperSingularCurvePairing1.CurveParams);
                                                    G2_1.SetCurveParams(SuperSingularCurvePairing1.CurveParams);
                                                    G2_1:=SuperSingularCurvePairing1.GetDefautGGenerator;
                                                  end;
                                                1:begin
                                                    SigAsPoint.SetCurveParams(MNTCurve1.CurveParams);
                                                    hashAspoint.SetCurveParams(MNTCurve1.CurveParams);
                                                    G2_2.SetCurveParams(MNTCurve1.CurveParams);
                                                    G2_2.SetToDefaultGenerator;
                                                  end;
                                                2:begin
                                                  SigAsPoint.SetCurveParams(BNCurvePairing1.CurveParams);
                                                  hashAspoint.SetCurveParams(BNCurvePairing1.CurveParams);
                                                  G2_2.SetCurveParams(BNCurvePairing1.CurveParams);
                                                  G2_2.SetToDefaultGenerator;
                                                  end;
                                                3:begin
                                                  SigAsPoint.SetCurveParams(BLS12CurvePairing1.CurveParams);
                                                  hashAspoint.SetCurveParams(BLS12CurvePairing1.CurveParams);
                                                  G2_2.SetCurveParams(BLS12CurvePairing1.CurveParams);
                                                  G2_2.SetToDefaultGenerator;
                                                  end;
                                                4:begin
                                                  SigAsPoint.SetCurveParams(BLS24CurvePairing1.CurveParams);
                                                  hashAspoint.SetCurveParams(BLS24CurvePairing1.CurveParams);
                                                  G2_4.SetCurveParams(BLS24CurvePairing1.CurveParams);
                                                  G2_4.SetToDefaultGenerator;
                                                  end;
                                                5:begin
                                                  SigAsPoint.SetCurveParams(BLS27Curve1.CurveParams);
                                                  hashAspoint.SetCurveParams(BLS27Curve1.CurveParams);
                                                  G2_9.SetCurveParams(BLS27Curve1.CurveParams);
                                                  G2_9.SetToDefaultGenerator;
                                                  end;
                                                6:begin
                                                  SigAsPoint.SetCurveParams(KSS16CurvePairing1.CurveParams);
                                                  hashAspoint.SetCurveParams(KSS16CurvePairing1.CurveParams);
                                                  G2_4.SetCurveParams(KSS16CurvePairing1.CurveParams);
                                                  G2_4.SetToDefaultGenerator;
                                                  end;
                                                7:begin
                                                  SigAsPoint.SetCurveParams(KSS18CurvePairing1.CurveParams);
                                                  hashAspoint.SetCurveParams(KSS18CurvePairing1.CurveParams);
                                                  G2_3.SetCurveParams(KSS18CurvePairing1.CurveParams);
                                                  G2_3.SetToDefaultGenerator;
                                                  end;
                                                8:begin
                                                  SigAsPoint.SetCurveParams(KSS36Curve1.CurveParams);
                                                  hashAspoint.SetCurveParams(KSS36Curve1.CurveParams);
                                                  G2_6.SetCurveParams(KSS36Curve1.CurveParams);
                                                  G2_6.SetToDefaultGenerator;
                                                  end;
                                                end;
                                                FileHash:=SHA256FileHash(filename);
                                                hashAspoint.SetAsTorsionFromHash(FileHash);
                                                FileSignature:=TFile.ReadAllBytes(Edit26.Text);
                                                SigAsPoint.DeCompressFromArray(FileSignature);
                                                case ComboBox10.ItemIndex of
                                                0:begin
                                                  e1_2:=SuperSingularCurvePairing1.Paire(SigAsPoint,G2_1);
                                                  e2_2:=SuperSingularCurvePairing1.Paire(hashAspoint,PubKeyInEFp);
                                                  Equals:=e1_12=e2_12;
                                                  end;
                                                1:begin
                                                  e1_6:=MNTCurve1.Paire(SigAsPoint,G2_3);
                                                  e2_6:=MNTCurve1.Paire(hashAspoint,PubKeyInEFp3);
                                                  Equals:=e1_6=e2_6;
                                                  end;
                                                2:begin
                                                  e1_12:=BNCurvePairing1.Paire(SigAsPoint,G2_2);
                                                  e2_12:=BNCurvePairing1.Paire(hashAspoint,PubKeyInEFp2);
                                                  Equals:=e1_12=e2_12;
                                                  end;
                                                3:begin
                                                  e1_12:=BLS12CurvePairing1.Paire(SigAsPoint,G2_2);
                                                  e2_12:=BLS12CurvePairing1.Paire(hashAspoint,PubKeyInEFp2);
                                                  Equals:=e1_12=e2_12;
                                                  end;
                                                4:begin
                                                  e1_24:=BLS24CurvePairing1.Paire(SigAsPoint,G2_4);
                                                  e2_24:=BLS24CurvePairing1.Paire(hashAspoint,PubKeyInEFp4);
                                                  Equals:=e1_12=e2_12;
                                                  end;
                                                5:begin
                                                  e1_27:=BLS27Curve1.Paire(SigAsPoint,G2_9);
                                                  e2_27:=BLS27Curve1.Paire(hashAspoint,PubKeyInEFp9);
                                                  Equals:=e1_27=e2_27;
                                                  end;
                                                6:begin
                                                  e1_16:=KSS16CurvePairing1.Paire(SigAsPoint,G2_4);
                                                  e2_16:=KSS16CurvePairing1.Paire(hashAspoint,PubKeyInEFp4);
                                                  Equals:=e1_16=e2_16;
                                                  end;
                                                7:begin
                                                  e1_18:=KSS18CurvePairing1.Paire(SigAsPoint,G2_3);
                                                  e2_18:=KSS18CurvePairing1.Paire(hashAspoint,PubKeyInEFp3);
                                                  Equals:=e1_18=e2_18;
                                                  end;
                                                8:begin
                                                  e1_36:=KSS36Curve1.Paire(SigAsPoint,G2_6);
                                                  e2_36:=KSS36Curve1.Paire(hashAspoint,PubKeyInEFp6);
                                                  Equals:=e1_36=e2_36;
                                                  end;
                                                end;
                                                if  Equals then MessageDlg('Signature Valide :'+ArrayToHex(FileSignature),mtInformation,[mbok],0)
                                                else Messagedlg('Signature Invalide :la signature ne correspond pas à la clé de vérification, ou bien au fichier signé',mterror,[mbok],0);
                                                end
                                else MessageDlg('Fichier de signature introuvable ',mtwarning,[mbok],0);
                                end
else MessageDlg('Fichier source à vérifié introuvable ',mtwarning,[mbok],0);
end;

procedure TForm2.Button4Click(Sender: TObject);
var MessageHash,FileSignature:TBytes;
    SigAsPoint:FpPoint;
begin
MessageHash:=SHA256StringHash(Memo1.Text);
case ComboBox10.ItemIndex of
0:SigAsPoint.SetCurveParams(SuperSingularCurvePairing1.CurveParams);
1:SigAsPoint.SetCurveParams(MNTCurve1.CurveParams);
2:SigAsPoint.SetCurveParams(BNCurvePairing1.CurveParams);
3:SigAsPoint.SetCurveParams(BLS12CurvePairing1.CurveParams);
4:SigAsPoint.SetCurveParams(BLS24CurvePairing1.CurveParams);
5:SigAsPoint.SetCurveParams(BLS27Curve1.CurveParams);
6:SigAsPoint.SetCurveParams(KSS16CurvePairing1.CurveParams);
7:SigAsPoint.SetCurveParams(KSS18CurvePairing1.CurveParams);
8:SigAsPoint.SetCurveParams(KSS36Curve1.CurveParams);
end;
SigAsPoint.SetAsTorsionFromHash(MessageHash);
SigAsPoint:=SecretKey*SigAsPoint;
FileSignature:=SigAsPoint.CompressToArray;
Edit1.Text:=ArrayToHex(FileSignature);
Label46.Caption:='('+inttostr((length(FileSignature)-1)*8)+' bit)';
MessageDlg('Signature Générée avec succés !',mtConfirmation,[mbok],0);
end;

procedure TForm2.Button5Click(Sender: TObject);
var FileSignature,Filehash:TBytes;
    SigAsPoint:FpPoint;
    hashAspoint:FpPoint;
    Equals:boolean;
begin
case ComboBox10.ItemIndex of
  0:begin
  SigAsPoint.SetCurveParams(SuperSingularCurvePairing1.CurveParams);
  hashAspoint.SetCurveParams(SuperSingularCurvePairing1.CurveParams);
  G2_1.SetCurveParams(SuperSingularCurvePairing1.CurveParams);
  G2_1:=SuperSingularCurvePairing1.GetDefautGGenerator;
  end;
  1:begin
  SigAsPoint.SetCurveParams(MNTCurve1.CurveParams);
  hashAspoint.SetCurveParams(MNTCurve1.CurveParams);
  G2_2.SetCurveParams(MNTCurve1.CurveParams);
  G2_2.SetToDefaultGenerator;
  end;
  2:begin
  SigAsPoint.SetCurveParams(BNCurvePairing1.CurveParams);
  hashAspoint.SetCurveParams(BNCurvePairing1.CurveParams);
  G2_2.SetCurveParams(BNCurvePairing1.CurveParams);
  G2_2.SetToDefaultGenerator;
  end;
  3:begin
  SigAsPoint.SetCurveParams(BLS12CurvePairing1.CurveParams);
  hashAspoint.SetCurveParams(BLS12CurvePairing1.CurveParams);
  G2_2.SetCurveParams(BLS12CurvePairing1.CurveParams);
  G2_2.SetToDefaultGenerator;
  end;
  4:begin
  SigAsPoint.SetCurveParams(BLS24CurvePairing1.CurveParams);
  hashAspoint.SetCurveParams(BLS24CurvePairing1.CurveParams);
  G2_4.SetCurveParams(BLS24CurvePairing1.CurveParams);
  G2_4.SetToDefaultGenerator;
  end;
  5:begin
  SigAsPoint.SetCurveParams(BLS27Curve1.CurveParams);
  hashAspoint.SetCurveParams(BLS27Curve1.CurveParams);
  G2_9.SetCurveParams(BLS27Curve1.CurveParams);
  G2_9.SetToDefaultGenerator;
  end;
  6:begin
  SigAsPoint.SetCurveParams(KSS16CurvePairing1.CurveParams);
  hashAspoint.SetCurveParams(KSS16CurvePairing1.CurveParams);
  G2_4.SetCurveParams(KSS16CurvePairing1.CurveParams);
  G2_4.SetToDefaultGenerator;
  end;
  7:begin
  SigAsPoint.SetCurveParams(KSS18CurvePairing1.CurveParams);
  hashAspoint.SetCurveParams(KSS18CurvePairing1.CurveParams);
  G2_3.SetCurveParams(KSS18CurvePairing1.CurveParams);
  G2_3.SetToDefaultGenerator;
  end;
  8:begin
  SigAsPoint.SetCurveParams(KSS36Curve1.CurveParams);
  hashAspoint.SetCurveParams(KSS36Curve1.CurveParams);
  G2_6.SetCurveParams(KSS36Curve1.CurveParams);
  G2_6.SetToDefaultGenerator;
  end;
  end;
  FileHash:=SHA256StringHash(Memo1.Text);
  hashAspoint.SetAsTorsionFromHash(FileHash);
  FileSignature:=TFile.ReadAllBytes(Edit26.Text);
  SigAsPoint.DeCompressFromArray(FileSignature);
  case ComboBox10.ItemIndex of
  0:begin
  e1_2:=SuperSingularCurvePairing1.Paire(SigAsPoint,G2_1);
  e2_2:=SuperSingularCurvePairing1.Paire(hashAspoint,PubKeyInEFp);
  Equals:=e1_12=e2_12;
  end;
  1:begin
  e1_6:=MNTCurve1.Paire(SigAsPoint,G2_3);
  e2_6:=MNTCurve1.Paire(hashAspoint,PubKeyInEFp3);
  Equals:=e1_6=e2_6;
  end;
  2:begin
  e1_12:=BNCurvePairing1.Paire(SigAsPoint,G2_2);
  e2_12:=BNCurvePairing1.Paire(hashAspoint,PubKeyInEFp2);
  Equals:=e1_12=e2_12;
  end;
  3:begin
  e1_12:=BLS12CurvePairing1.Paire(SigAsPoint,G2_2);
  e2_12:=BLS12CurvePairing1.Paire(hashAspoint,PubKeyInEFp2);
  Equals:=e1_12=e2_12;
  end;
  4:begin
  e1_24:=BLS24CurvePairing1.Paire(SigAsPoint,G2_4);
  e2_24:=BLS24CurvePairing1.Paire(hashAspoint,PubKeyInEFp4);
  Equals:=e1_24=e2_24;
  end;
  5:begin
  e1_27:=BLS27Curve1.Paire(SigAsPoint,G2_9);
  e2_27:=BLS27Curve1.Paire(hashAspoint,PubKeyInEFp9);
  Equals:=e1_27=e2_27;
  end;
  6:begin
  e1_16:=KSS16CurvePairing1.Paire(SigAsPoint,G2_4);
  e2_16:=KSS16CurvePairing1.Paire(hashAspoint,PubKeyInEFp4);
  Equals:=e1_16=e2_16;
  end;
  7:begin
  e1_18:=KSS18CurvePairing1.Paire(SigAsPoint,G2_3);
  e2_18:=KSS18CurvePairing1.Paire(hashAspoint,PubKeyInEFp3);
  Equals:=e1_18=e2_18;
  end;
  8:begin
  e1_36:=KSS36Curve1.Paire(SigAsPoint,G2_6);
  e2_36:=KSS36Curve1.Paire(hashAspoint,PubKeyInEFp6);
  Equals:=e1_36=e2_36;
  end;
  end;
  if  Equals then MessageDlg('Signature Valide :'+ArrayToHex(FileSignature),mtInformation,[mbok],0)
  else Messagedlg('Signature Invalide :la signature ne correspond pas à la clé de vérification, ou bien au fichier signé',mterror,[mbok],0);
end;

procedure TForm2.ComboBox10Change(Sender: TObject);
var i:integer;
begin
ComboBox1.Items.Clear;
ComboBox11.Items.Clear;
case
  ComboBox10.ItemIndex of
  0:for i:=0 to length(SSParamsList)-1 do ComboBox1.Items.add(SSparamsList[i]);
  1:for i:=0 to length(MNTParamsList)-1 do ComboBox1.Items.add(MNTparamsList[i]);
  2:for i:=0 to length(BNParamsList)-1 do ComboBox1.Items.add(bnparamsList[i]);
  3:for i:=0 to length(BLS12ParamsList)-1 do ComboBox1.Items.add(BLS12paramsList[i]);
  4:for i:=0 to length(BLS24ParamsList)-1 do ComboBox1.Items.add(BLS24paramsList[i]);
  5:for i:=0 to length(BLS27ParamsList)-1 do ComboBox1.Items.add(BLS27paramsList[i]);
  6:for i:=0 to length(KSS16ParamsList)-1 do ComboBox1.Items.add(KSS16paramsList[i]);
  7:for i:=0 to length(KSS18ParamsList)-1 do ComboBox1.Items.add(KSS18paramsList[i]);
  8:for i:=0 to length(KSS36ParamsList)-1 do ComboBox1.Items.add(KSS36paramsList[i]);
end;
case
  ComboBox10.ItemIndex of
  0:for i:=0 to length(SSImplementedPairingAlgos)-1 do ComboBox11.Items.add(SSImplementedPairingAlgos[i]);
  1:for i:=0 to length(MNTImplementedPairingAlgos)-1 do ComboBox11.Items.add(MNTImplementedPairingAlgos[i]);
  2:for i:=0 to length(BNImplementedPairingAlgos)-1 do ComboBox11.Items.add(BNImplementedPairingAlgos[i]);
  3:for i:=0 to length(BLS12ImplementedPairingAlgos)-1 do ComboBox11.Items.add(BLS12ImplementedPairingAlgos[i]);
  4:for i:=0 to length(BLS24ImplementedPairingAlgos)-1 do ComboBox11.Items.add(BLS24ImplementedPairingAlgos[i]);
  5:for i:=0 to length(BLS27ImplementedPairingAlgos)-1 do ComboBox11.Items.add(BLS27ImplementedPairingAlgos[i]);
  6:for i:=0 to length(KSS16ImplementedPairingAlgos)-1 do ComboBox11.Items.add(KSS16ImplementedPairingAlgos[i]);
  7:for i:=0 to length(KSS18ImplementedPairingAlgos)-1 do ComboBox11.Items.add(KSS18ImplementedPairingAlgos[i]);
  8:for i:=0 to length(KSS36ImplementedPairingAlgos)-1 do ComboBox11.Items.add(KSS36ImplementedPairingAlgos[i]);
end;
case
  ComboBox10.ItemIndex of
  0:begin Label32.Caption:='Public Key (a point on E(Fp2))';Label32.Caption:='E(Fp2) Default Generator ';end;
  1:begin Label32.Caption:='Public Key (a point on E(Fp3))';Label32.Caption:='E(Fp3) Default Generator ';end;
  2:begin Label32.Caption:='Public Key (a point on E(Fp2))';Label32.Caption:='E(Fp2) Default Generator ';end;
  3:begin Label32.Caption:='Public Key (a point on E(Fp2))';Label32.Caption:='E(Fp2) Default Generator ';end;
  4:begin Label32.Caption:='Public Key (a point on E(Fp4))';Label32.Caption:='E(Fp4) Default Generator ';end;
  5:begin Label32.Caption:='Public Key (a point on E(Fp9))';Label32.Caption:='E(Fp9) Default Generator ';end;
  6:begin Label32.Caption:='Public Key (a point on E(Fp4))';Label32.Caption:='E(Fp4) Default Generator ';end;
  7:begin Label32.Caption:='Public Key (a point on E(Fp3))';Label32.Caption:='E(Fp3) Default Generator ';end;
  8:begin Label32.Caption:='Public Key (a point on E(Fp6))';Label32.Caption:='E(Fp6) Default Generator ';end;
end;
ComboBox1.ItemIndex:=0;
ComboBox11.ItemIndex:=0;
end;

procedure TForm2.ComboBox11Change(Sender: TObject);
begin
case ComboBox10.ItemIndex of
0:SuperSingularCurvePairing1.PairingAlgorithm:=SSCurvesPairingAlgos(ComboBox11.ItemIndex);
2:MNTCurve1.PairingAlgorithm:=MNTCurvesPairingAlgos(ComboBox11.ItemIndex);
3:BNCurvePairing1.PairingAlgorithm:=BnCurvesPairingAlgos(ComboBox11.ItemIndex);
4:BLS12CurvePairing1.PairingAlgorithm:=BLS12CurvesPairingAlgos(ComboBox11.ItemIndex);
5:BLS24CurvePairing1.PairingAlgorithm:=BLS24CurvesPairingAlgos(ComboBox11.ItemIndex);
6:BLS27Curve1.PairingAlgorithm:=BLS27CurvesPairingAlgos(ComboBox11.ItemIndex);
7:KSS16CurvePairing1.PairingAlgorithm:=KSS16CurvesPairingAlgos(ComboBox11.ItemIndex);
8:KSS18CurvePairing1.PairingAlgorithm:=KSS18CurvesPairingAlgos(ComboBox11.ItemIndex);
9:KSS36Curve1.PairingAlgorithm:=KSS36CurvesPairingAlgos(ComboBox11.ItemIndex);
end;
end;

procedure TForm2.ComboBox1Change(Sender: TObject);
begin
case ComboBox10.ItemIndex of
0:SuperSingularCurvePairing1.Parametres:=StandardSSCurves(ComboBox1.ItemIndex);
2:MNTCurve1.Parametres:=StandardMNTCurves(ComboBox1.ItemIndex);
3:BNCurvePairing1.Parametres:=StandardBNCurves(ComboBox1.ItemIndex);
4:BLS12CurvePairing1.Parametres:=StandardBLS12Curves(ComboBox1.ItemIndex);
5:BLS24CurvePairing1.Parametres:=StandardBLS24Curves(ComboBox1.ItemIndex);
7:KSS16CurvePairing1.Parametres:=StandardKSS16Curves(ComboBox1.ItemIndex);
8:KSS18CurvePairing1.Parametres:=StandardKSS18Curves(ComboBox1.ItemIndex);
9:KSS36Curve1.Parametres:=StandardKSS36Curves(ComboBox1.ItemIndex);
end;
end;

procedure TForm2.FormShow(Sender: TObject);
begin
ComboBox10Change(Self);
ComboBox1.ItemIndex:=7;
ComboBox11.ItemIndex:=1;
end;

procedure TForm2.SpeedButton1Click(Sender: TObject);
begin
if OpenDialog1.Execute then Edit24.Text:=OpenDialog1.FileName;
end;

procedure TForm2.SpeedButton2Click(Sender: TObject);
begin
if SaveDialog1.Execute then Edit27.Text:=SaveDialog1.FileName;
end;

procedure TForm2.SpeedButton3Click(Sender: TObject);
begin
if OpenDialog1.Execute then edit25.Text:=OpenDialog1.FileName;
end;

procedure TForm2.SpeedButton4Click(Sender: TObject);
begin
if OpenDialog1.Execute then edit26.Text:=OpenDialog1.FileName;
end;

procedure TForm2.SpeedButton5Click(Sender: TObject);
var filekey:string;
    list:TStringList;
    g:Fp2Point;
    tmp:TBytes;
begin
list:=TStringList.Create;
if OpenDialog1.Execute then begin
                            list.LoadFromFile(OpenDialog1.FileName);
                            try
                            ComboBox10.ItemIndex:=StrToInt(List.Strings[0]);
                            ComboBox1.ItemIndex:=StrToInt(List.Strings[1]);
                            ComboBox11.ItemIndex:=StrToInt(List.Strings[2]);
                            ComboBox10Change(Self);
                            ComboBox1Change(Self);
                            ComboBox11Change(Self);
                            SecretKey:=List.Strings[3];
                            tmp:=HexToArray(List.Strings[4]);
                            Edit31.Text:=List.Strings[4];
                            case ComboBox10.ItemIndex of
                            0:begin     // Super Singular
                                PubKeyInEFp.SetCurveParams(SuperSingularCurvePairing1.CurveParams);
                                PubKeyInEFp:=SuperSingularCurvePairing1.GetDefautGGenerator;
                                Edit30.Text:=PubKeyInEFp.ToHexString;
                                PubKeyInEFp.DeCompressFromArray(tmp);
                                Edit23.Text:=PubKeyInEFp.toHexString;
                              end;
                            1:begin    // MNT
                                PubKeyInEFp3.SetCurveParams(MNTCurve1.CurveParams);
                                PubKeyInEFp3.SetToDefaultGenerator;
                                Edit30.Text:=PubKeyInEFp3.ToHexString;
                                PubKeyInEFp3.DeCompressFromArray(tmp);
                                Edit23.Text:=PubKeyInEFp.toHexString;
                              end;
                            2:begin     //BN
                                PubKeyInEFp2.SetCurveParams(BNCurvePairing1.CurveParams);
                                PubKeyInEFp2.SetToDefaultGenerator;
                                Edit30.Text:=PubKeyInEFp2.ToHexString;
                                PubKeyInEFp2.DeCompressFromArray(tmp);
                                Edit23.Text:=PubKeyInEFp2.toHexString;
                              end;
                            3:begin      //BLS12
                                PubKeyInEFp2.SetCurveParams(BLS12CurvePairing1.CurveParams);
                                PubKeyInEFp2.SetToDefaultGenerator;
                                Edit30.Text:=PubKeyInEFp2.ToHexString;
                                PubKeyInEFp2.DeCompressFromArray(tmp);
                                Edit23.Text:=PubKeyInEFp2.toHexString;
                              end;
                            4:begin       //BLS24
                                PubKeyInEFp4.SetCurveParams(BLS24CurvePairing1.CurveParams);
                                PubKeyInEFp4.SetToDefaultGenerator;
                                Edit30.Text:=PubKeyInEFp4.ToHexString;
                                PubKeyInEFp4.DeCompressFromArray(tmp);
                                Edit23.Text:=PubKeyInEFp4.toHexString;
                              end;
                            5:begin       //BLS27
                                PubKeyInEFp9.SetCurveParams(BLS27Curve1.CurveParams);
                                PubKeyInEFp9.SetToDefaultGenerator;
                                Edit30.Text:=PubKeyInEFp9.ToHexString;
                                PubKeyInEFp9.DeCompressFromArray(tmp);
                                Edit23.Text:=PubKeyInEFp9.toHexString;
                              end;
                            6:begin       //KSS16
                                PubKeyInEFp4.SetCurveParams(KSS16CurvePairing1.CurveParams);
                                PubKeyInEFp4.SetToDefaultGenerator;
                                Edit30.Text:=PubKeyInEFp4.ToHexString;
                                PubKeyInEFp4.DeCompressFromArray(tmp);
                                Edit23.Text:=PubKeyInEFp4.toHexString;
                              end;
                            7:begin      //KSS18
                                PubKeyInEFp3.SetCurveParams(KSS18CurvePairing1.CurveParams);
                                PubKeyInEFp3.SetToDefaultGenerator;
                                Edit30.Text:=PubKeyInEFp3.ToHexString;
                                PubKeyInEFp3.DeCompressFromArray(tmp);
                                Edit23.Text:=PubKeyInEFp3.toHexString;

                            end;
                            8:begin       //KSS36
                                PubKeyInEFp6.SetCurveParams(KSS36Curve1.CurveParams);
                                PubKeyInEFp6.SetToDefaultGenerator;
                                Edit30.Text:=PubKeyInEFp6.ToHexString;
                                PubKeyInEFp6.DeCompressFromArray(tmp);
                                Edit23.Text:=PubKeyInEFp6.toHexString;
                              end;
                            end;
                            label31.Caption:='Secrete Key (an element from FP) : on '+inttostr(SecretKey.BitLength)+' bit';
                            Edit22.Text:=SecretKey.ToHexString;
                            Button2.Enabled:=true;
                            Button3.Enabled:=true;
                            Button4.Enabled:=true;
                            Button5.Enabled:=true;
                            SpeedButton5.Enabled:=true;
                            SpeedButton6.Enabled:=true;
                            except
                            MessageDlg('Fichier contenant la paire des clés invalide ',mterror,[mbok],0);
                            end;
                            end;
end;

procedure TForm2.SpeedButton6Click(Sender: TObject);
var list:TStringList;
    filekey:string;
begin
list:=TStringList.Create;
list.Add(inttostr(ComboBox10.ItemIndex));
list.Add(inttostr(ComboBox1.ItemIndex));
list.Add(inttostr(ComboBox11.ItemIndex));
list.Add(edit22.Text);
List.Add(edit31.Text);
if SaveDialog1.Execute then begin
                            filekey:=SaveDialog1.FileName;
                            if Pos('.key',LowerCase(filekey))=0 then filekey:=filekey+'.key';
                            list.SaveToFile(filekey);
                            end;
end;

end.
