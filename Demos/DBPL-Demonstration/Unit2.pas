unit Unit2;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, BNCurves, KSS36Curves, BLS27Curves,
  KSS18Curves, SuperSingularCurves, KSS16Curves, BLS24Curves, BLS12Curves,
  Vcl.ComCtrls, Vcl.ExtCtrls, GeneralTypes, MNTCurves, System.ImageList,
  Vcl.ImgList, Vcl.Menus, Vcl.StdCtrls, TreeView1, Vcl.ToolWin;

type
  TForm2 = class(TForm)
    Panel1: TPanel;
    TreeView1: TTreeView;
    BNCurvePairing1: TBNCurvePairing;
    ImageList1: TImageList;
    MainMenu1: TMainMenu;
    About1: TMenuItem;
    Splitter1: TSplitter;
    GroupBox1: TGroupBox;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Splitter2: TSplitter;
    ComboBox1: TComboBox;
    ComboBox2: TComboBox;
    ComboBox3: TComboBox;
    ComboBox4: TComboBox;
    Button1: TButton;
    BLS24CurvePairing1: TBLS24CurvePairing;
    BLS12CurvePairing1: TBLS12CurvePairing;
    KSS16CurvePairing1: TKSS16CurvePairing;
    MNTCurve1: TMNTCurve;
    SuperSingularCurvePairing1: TSuperSingularCurvePairing;
    KSS18CurvePairing1: TKSS18CurvePairing;
    BLS27Curve1: TBLS27Curve;
    KSS36Curve1: TKSS36Curve;
    Panel2: TPanel;
    TreeView11: TTreeView1;
    Memo1: TRichEdit;
    procedure FormShow(Sender: TObject);
    procedure TreeView1Click(Sender: TObject);
    procedure Button1Click(Sender: TObject);
    procedure Memo1Change(Sender: TObject);
    procedure TreeView11CustomDrawItem(Sender: TCustomTreeView; Node: TTreeNode;
      State: TCustomDrawState; var DefaultDraw: Boolean);
    procedure About1Click(Sender: TObject);
  private
    { Déclarations privées }
  public
  procedure Print(s1,s2: string;color:integer);
    { Déclarations publiques }
  end;

var
  Form2: TForm2;
  ActiveCurve:string;
  SelectedCurve:integer;
 /// Declare pairings variables here to avoid stack overflow ....
  var
  Pg1Bls12: G1BLS12;
  Pg2Bls12: G2BLS12;
  eBls12_1,eBls12_2: GTBLS12;

  Pg1Bls24: G1BLS24;
  Pg2Bls24: G2BLS24;
  eBls24_1,eBls24_2: GTBLS24;

  Pg1Bls27: G1BLS27;
  Pg2Bls27: G2BLS27;
  eBls27_1,eBls27_2: GTBLS27;

  Pg1Bn: G1Bn;
  Pg2Bn: G2Bn;
  eBn_1,eBn_2: GTBn;

  Pg1KSS16: G1Kss16;
  Pg2Kss16: G2Kss16;
  eKss16_1,eKss16_2: GTKss16;

  Pg1Kss18: G1Kss18;
  Pg2Kss18: G2Kss18;
  eKss18_1,eKss18_2: GTKss18;

  Pg1Kss36: G1Kss36;
  Pg2Kss36: G2Kss36;
  eKss36_1,eKss36_2: GTKss36;

  Pg1Mnt: G1Mnt;
  Pg2Mnt: G2Mnt;
  eMnt_1,eMnt_2: GTMnt;

  Pg1SS,Pg2SS: GSS;
  eSS_1,eSS_2:  GTSS;

  tim:TTime;
  h,m,s,ms:Word;
  time:String;
  a,b,c,ab:LargeInt;

implementation

{$R *.dfm}

uses Unit1;

procedure SSParms;
begin
Form2.ComboBox1.Items.Clear;
Form2.ComboBox1.Items.Add('Affine');
Form2.ComboBox1.Items.Add('Projective');
Form2.ComboBox1.Items.Add('Jacobian');
Form2.ComboBox1.ItemIndex:=1;
Form2.ComboBox2.Clear;
Form2.ComboBox2.Items.Add('Normal');
Form2.ComboBox2.ItemIndex:=0;
Form2.ComboBox3.Clear;
Form2.ComboBox3.Items.Add('Binary');
Form2.ComboBox3.Items.Add('Negative');
Form2.ComboBox3.ItemIndex:=1;
Form2.ComboBox4.Items.Clear;
Form2.ComboBox4.Items.Add('Tate');
Form2.ComboBox4.ItemIndex:=0;
end;

procedure BNParms;
begin
Form2.ComboBox1.Items.Clear;
Form2.ComboBox1.Items.Add('Affine');
Form2.ComboBox1.Items.Add('Projective');
Form2.ComboBox1.Items.Add('Jacobian');
Form2.ComboBox1.ItemIndex:=1;
Form2.ComboBox2.Clear;
Form2.ComboBox2.Items.Add('Normal');
Form2.ComboBox2.Items.Add('Karabina');
Form2.ComboBox2.ItemIndex:=1;
Form2.ComboBox2.Items.Add('Beuchat');
Form2.ComboBox3.Clear;
Form2.ComboBox3.Items.Add('Binary');
Form2.ComboBox3.Items.Add('Negative');
Form2.ComboBox3.ItemIndex:=1;
Form2.ComboBox4.Items.Clear;
Form2.ComboBox4.Items.Add('Optimal Ate');
Form2.ComboBox4.Items.Add('R-Ate');
Form2.ComboBox4.Items.Add('Tate');
Form2.ComboBox4.Items.Add('Eta');
Form2.ComboBox4.ItemIndex:=0;
end;
procedure BLS12Parms;
begin
Form2.ComboBox1.Items.Clear;
Form2.ComboBox1.Items.Add('Affine');
Form2.ComboBox1.Items.Add('Projective');
Form2.ComboBox1.Items.Add('Jacobian');
Form2.ComboBox1.ItemIndex:=1;
Form2.ComboBox2.Clear;
Form2.ComboBox2.Items.Add('Normal');
Form2.ComboBox2.Items.Add('Karabina');
Form2.ComboBox2.Items.Add('Beuchat');
Form2.ComboBox2.ItemIndex:=1;
Form2.ComboBox3.Clear;
Form2.ComboBox3.Items.Add('Binary');
Form2.ComboBox3.Items.Add('Negative');
Form2.ComboBox3.ItemIndex:=1;
Form2.ComboBox4.Items.Clear;
Form2.ComboBox4.Items.Add('Optimal Ate');
Form2.ComboBox4.ItemIndex:=0;
end;

procedure BLS24Parms;
begin
Form2.ComboBox1.Items.Clear;
Form2.ComboBox1.Items.Add('Affine');
Form2.ComboBox1.Items.Add('Projective');
Form2.ComboBox1.ItemIndex:=1;
Form2.ComboBox2.Clear;
Form2.ComboBox2.Items.Add('Normal');
Form2.ComboBox2.Items.Add('Karabina');
Form2.ComboBox2.ItemIndex:=1;
Form2.ComboBox3.Clear;
Form2.ComboBox3.Items.Add('Binary');
Form2.ComboBox3.Items.Add('Negative');
Form2.ComboBox3.ItemIndex:=1;
Form2.ComboBox4.Items.Clear;
Form2.ComboBox4.Items.Add('Optimal Ate');
Form2.ComboBox4.ItemIndex:=0;
end;

procedure BLS27Parms;
begin
Form2.ComboBox1.Items.Clear;
Form2.ComboBox1.Items.Add('Affine');
Form2.ComboBox1.ItemIndex:=0;
Form2.ComboBox2.Clear;
Form2.ComboBox2.Items.Add('Normal');
Form2.ComboBox2.ItemIndex:=1;
Form2.ComboBox3.Clear;
Form2.ComboBox3.Items.Add('Binary');
Form2.ComboBox3.Items.Add('Negative');
Form2.ComboBox3.ItemIndex:=2;
Form2.ComboBox4.Items.Clear;
Form2.ComboBox4.Items.Add('Optimal Ate');
Form2.ComboBox4.ItemIndex:=0;
end;

procedure KSS16Parms;
begin
Form2.ComboBox1.Items.Clear;
Form2.ComboBox1.Items.Add('Affine');
Form2.ComboBox1.Items.Add('Projective');
Form2.ComboBox1.ItemIndex:=1;
Form2.ComboBox2.Clear;
Form2.ComboBox2.Items.Add('Normal');
Form2.ComboBox2.ItemIndex:=0;
Form2.ComboBox3.Clear;
Form2.ComboBox3.Items.Add('Binary');
Form2.ComboBox3.Items.Add('Negative');
Form2.ComboBox3.ItemIndex:=1;
Form2.ComboBox4.Items.Clear;
Form2.ComboBox4.Items.Add('Optimal Ate');
Form2.ComboBox4.ItemIndex:=0;
end;

procedure kss18params;
begin
Form2.ComboBox1.Items.Clear;
Form2.ComboBox1.Items.Add('Affine');
Form2.ComboBox1.Items.Add('Projective');
Form2.ComboBox1.Items.Add('Jacobian');
Form2.ComboBox1.ItemIndex:=1;
Form2.ComboBox2.Clear;
Form2.ComboBox2.Items.Add('Normal');
Form2.ComboBox2.Items.Add('Karabina');
Form2.ComboBox2.Items.Add('Beuchat');
Form2.ComboBox2.ItemIndex:=1;
Form2.ComboBox3.Clear;
Form2.ComboBox3.Items.Add('Binary');
Form2.ComboBox3.Items.Add('Negative');
Form2.ComboBox3.ItemIndex:=1;
Form2.ComboBox4.Items.Clear;
Form2.ComboBox4.Items.Add('Optimal Ate');
Form2.ComboBox4.ItemIndex:=0;
end;

procedure kss36params;
begin
Form2.ComboBox1.Items.Clear;
Form2.ComboBox1.Items.Add('Affine');
Form2.ComboBox1.Items.Add('Projective');
Form2.ComboBox1.ItemIndex:=1;
Form2.ComboBox2.Clear;
Form2.ComboBox2.Items.Add('Normal');
Form2.ComboBox2.ItemIndex:=0;
Form2.ComboBox3.Clear;
Form2.ComboBox3.Items.Add('Binary');
Form2.ComboBox3.Items.Add('Negative');
Form2.ComboBox3.ItemIndex:=1;
Form2.ComboBox4.Items.Clear;
Form2.ComboBox4.Items.Add('Optimal Ate');
Form2.ComboBox4.ItemIndex:=0;
end;

procedure Mntparams;
begin
Form2.ComboBox1.Items.Clear;
Form2.ComboBox1.Items.Add('Affine');
Form2.ComboBox1.Items.Add('Projective');
Form2.ComboBox1.ItemIndex:=0;
Form2.ComboBox2.Clear;
Form2.ComboBox2.Items.Add('Normal');
Form2.ComboBox2.ItemIndex:=0;
Form2.ComboBox3.Clear;
Form2.ComboBox3.Items.Add('Binary');
Form2.ComboBox3.Items.Add('Negative');
Form2.ComboBox3.ItemIndex:=0;
Form2.ComboBox4.Items.Clear;
Form2.ComboBox4.Items.Add('Tate');
Form2.ComboBox4.Items.Add('Optimal Ate');
Form2.ComboBox4.ItemIndex:=0;
end;

procedure TForm2.Print(s1,s2: string;color:integer);
begin
Memo1.SelStart:=length(Memo1.Text);
Memo1.SelLength:=length(s1);
Memo1.SelAttributes.Style:=Memo1.SelAttributes.Style+[fsbold];
Memo1.Lines.Add(s1);
Memo1.SelAttributes.Style:=Memo1.SelAttributes.Style-[fsbold];
case color of
1:Memo1.SelAttributes.Color:=clGreen;
2:begin
  Memo1.SelAttributes.Color:=clMaroon;
  Memo1.SelAttributes.Style:=Memo1.SelAttributes.Style+[fsbold];
  end;
else Memo1.SelAttributes.Color:=clBlack;
end;
Memo1.Lines.Add('   '+s2);
Memo1.Lines.Add('   ');
Application.ProcessMessages;
end;

procedure TForm2.About1Click(Sender: TObject);
begin
form1.showmodal;
end;

procedure TForm2.Button1Click(Sender: TObject);
begin
if SelectedCurve=6 then begin
                        case combobox1.ItemIndex of
                        0:BNCurvePairing1.CoordinatesSystem:=csAffine;
                        1:BNCurvePairing1.CoordinatesSystem:=csProjective;
                        2:BNCurvePairing1.CoordinatesSystem:=csJacobian;
                        end;
                        case combobox2.ItemIndex of
                        0:BNCurvePairing1.Fp12PoweringMode:=pmNormal;
                        1:BNCurvePairing1.Fp12PoweringMode:=pmKarbina;
                        2:BNCurvePairing1.Fp12PoweringMode:=pmBeuchat;
                        end;
                        case combobox4.ItemIndex of
                        0:BNCurvePairing1.PairingAlgorithm:=bpOptAte;
                        1:BNCurvePairing1.PairingAlgorithm:=bpR_Ate;
                        2:BNCurvePairing1.PairingAlgorithm:=bpEta;
                        3:BNCurvePairing1.PairingAlgorithm:=bpTate;
                        end;
                        case combobox3.ItemIndex of
                        0:BNCurvePairing1.LoopMode:=lpmBinary;
                        1:BNCurvePairing1.LoopMode:=lpmNaf;
                        2:BNCurvePairing1.LoopMode:=lpmAuto;
                        end;
                        Print('BN Curves Bilinearity Test :('+TreeView1.Selected.Text+' Curve.)' ,'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'----------------------------------------------------------------------------------------------------------------------------------',0);
                        Pg1Bn:=BNCurvePairing1.HashToG1Point('kamel');
                        Pg2Bn:=BNCurvePairing1.HashToG2Point('kamel');
                        Print('First Point P (from E(Fp)):',Pg1Bn.toHexString,0);
                        Print('Second Point Q (from E(Fp2)):',Pg2BN.ToHexString,0);
                        tim:=now;
                        ebn_1:=BNCurvePairing1.Paire(Pg1BN,Pg2Bn);
                        tim:=now-tim;
                        DecodeTime(tim,h,m,s,ms);
                        time:='(in '+inttostr(m)+':'+inttostr(s)+':'+inttostr(ms)+'ms)';
                        Print('Pairing e(P,Q): '+time,ebn_1.toHexString,0);
                        a.SetToRandom(BNCurvePairing1.P);
                        a:=10;
                        Print('First Scalar Exponent a (From Fp):',a.toHexString,0);
                        b.SetToRandom(BNCurvePairing1.P);
                        b:=20;
                        Print('Second Scalar Exponent b (From Fp):',b.toHexString,0);
                        ab:=a*b;
                        Print('Product a*b:',ab.toHexString,0);
                        ebn_1:=ebn_1.Pow(ab,pmKarbina);
                        Print('Pairing''s exponent e(P,Q)^(a*b):',ebn_1.toHexString,1);
                        Pg1bn:=a*Pg1bn;
                        Pg2Bn:=b*Pg2bn;
                        Print('Point a*P :',pg1bn.toHexString,7);
                        Print('Point b*Q :',pg2bn.toHexString,8);
                        ebn_2:=BNCurvePairing1.Paire(pg1bn,pg2bn);
                        Print('Pairing e(aP,bQ): ',ebn_2.toHexString,1);
                        if eBn_2=ebn_1 then Print('','Bilinearity Verified e(P,Q)^(a*b)= e(aP,bQ)' ,2)
                        else Print('','Bilinearity not Verified, somthing wrong' ,2);
                        Print('','-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'-------------------------------------------------------------------------------------------------------',0);
                        end
else if SelectedCurve=0 then begin
                        Print('BLS12 Curves Bilinearity Test :('+TreeView1.Selected.Text+' Curve.)' ,'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'----------------------------------------------------------------------------------------------------------------------------------',0);
                        Pg1Bls12:=BLS12CurvePairing1.HashToG1Point('kamel');
                        Pg2Bls12:=BLS12CurvePairing1.HashToG2Point('kamel');
                        Print('First Point P (from E(Fp)):',Pg1Bls12.toHexString,0);
                        Print('Second Point Q (from E(Fp2)):',Pg2Bls12.ToHexString,0);
                        tim:=now;
                        eBls12_1:=BLS12CurvePairing1.Paire(Pg1Bls12,Pg2Bls12);
                        tim:=now-tim;
                        DecodeTime(tim,h,m,s,ms);
                        time:='(in '+inttostr(m)+':'+inttostr(s)+':'+inttostr(ms)+'ms)';
                        Print('Pairing e(P,Q): '+time,eBls12_1.toHexString,0);
                        a.SetToRandom(BLS12CurvePairing1.P);
                        Print('First Scalar Exponent a (From Fp):',a.toHexString,0);
                        b.SetToRandom(BLS12CurvePairing1.P);
                        Print('Second Scalar Exponent b (From Fp):',b.toHexString,0);
                        ab:=a*b;
                        Print('Product a*b:',ab.toHexString,0);
                        eBls12_1:=eBls12_1.Pow(ab,pmKarbina);
                        Print('Pairing''s exponent e(P,Q)^(a*b):',ebls12_1.toHexString,1);
                        Pg1Bls12:=a*Pg1Bls12;
                        Pg2Bls12:=b*Pg2Bls12;
                        Print('Point a*P :',Pg1Bls12.toHexString,7);
                        Print('Point b*Q :',Pg2Bls12.toHexString,8);
                        ebls12_2:=BLs12CurvePairing1.Paire(Pg1Bls12,Pg2Bls12);
                        Print('Pairing e(aP,bQ): ',eBls12_2.toHexString,1);
                        if eBls12_1=eBls12_2 then Print('','Bilinearity Verified e(P,Q)^(a*b)= e(aP,bQ)' ,2)
                        else Print('','Bilinearity not Verified, somthing wrong' ,2);
                        Print('','-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'-------------------------------------------------------------------------------------------------------',0);
                        end
else if SelectedCurve=1 then begin
                        Print('Bls24 Curves Bilinearity Test :('+TreeView1.Selected.Text+' Curve.)' ,'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'----------------------------------------------------------------------------------------------------------------------------------',0);
                        Pg1Bls24:=Bls24CurvePairing1.HashToG1Point('kamel');
                        Pg2Bls24:=Bls24CurvePairing1.HashToG2Point('kamel');
                        Print('First Point P (from E(Fp)):',Pg1Bls24.toHexString,0);
                        Print('Second Point Q (from E(Fp2)):',Pg2Bls24.ToHexString,0);
                        tim:=now;
                        eBls24_1:=Bls24CurvePairing1.Paire(Pg1Bls24,Pg2Bls24);
                        tim:=now-tim;
                        DecodeTime(tim,h,m,s,ms);
                        time:='(in '+inttostr(m)+':'+inttostr(s)+':'+inttostr(ms)+'ms)';
                        Print('Pairing e(P,Q): '+time,eBls24_1.toHexString,0);
                        a.SetToRandom(Bls24CurvePairing1.P);
                        Print('First Scalar Exponent a (From Fp):',a.toHexString,0);
                        b.SetToRandom(BNCurvePairing1.P);
                        Print('Second Scalar Exponent b (From Fp):',b.toHexString,0);
                        ab:=a*b;
                        Print('Product a*b:',ab.toHexString,0);
                        eBls24_1:=eBls24_1.Pow(ab,pmKarbina);
                        Print('Pairing''s exponent e(P,Q)^(a*b):',eBls24_1.toHexString,1);
                        Pg1Bls24:=a*Pg1Bls24;
                        Pg2Bls24:=b*Pg2Bls24;
                        Print('Point a*P :',Pg1Bls24.toHexString,7);
                        Print('Point b*Q :',Pg2Bls24.toHexString,8);
                        eBls24_2:=Bls24CurvePairing1.Paire(Pg1Bls24,Pg2Bls24);
                        Print('Pairing e(aP,bQ): ',eBls24_2.toHexString,1);
                        if eBls24_1=eBls24_2 then Print('','Bilinearity Verified e(P,Q)^(a*b)= e(aP,bQ)' ,2)
                        else Print('','Bilinearity not Verified, somthing wrong' ,2);
                        Print('','-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'-------------------------------------------------------------------------------------------------------',0);
                        end
else if SelectedCurve=2 then begin
                        Print('Bls27 Curves Bilinearity Test :('+TreeView1.Selected.Text+' Curve.)' ,'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'----------------------------------------------------------------------------------------------------------------------------------',0);
                        Pg1Bls27:=Bls27Curve1.HashToG1Point('kamel');
                        Pg2Bls27:=Bls27Curve1.HashToG2Point('kamel');
                        Print('First Point P (from E(Fp)):',Pg1Bls27.toHexString,0);
                        Print('Second Point Q (from E(Fp2)):',Pg2Bls27.ToHexString,0);
                        tim:=now;
                        eBls27_1:=BLS27Curve1.Paire(Pg1Bls27,Pg2Bls27);
                        tim:=now-tim;
                        DecodeTime(tim,h,m,s,ms);
                        time:='(in '+inttostr(m)+':'+inttostr(s)+':'+inttostr(ms)+'ms)';
                        Print('Pairing e(P,Q): '+time,eBls27_1.toHexString,0);
                        a.SetToRandom(Bls27Curve1.P);
                        Print('First Scalar Exponent a (From Fp):',a.toHexString,0);
                        b.SetToRandom(Bls27Curve1.P);
                        Print('Second Scalar Exponent b (From Fp):',b.toHexString,0);
                        ab:=a*b;
                        Print('Product a*b:',ab.toHexString,0);
                        eBls27_1:=eBls27_1.Pow(ab);
                        Print('Pairing''s exponent e(P,Q)^(a*b):',eBls27_1.toHexString,1);
                        Pg1Bls27:=a*Pg1Bls27;
                        Pg2Bls27:=b*Pg2Bls27;
                        Print('Point a*P :',Pg1Bls27.toHexString,7);
                        Print('Point b*Q :',Pg2Bls27.toHexString,8);
                        eBls27_2:=Bls27Curve1.Paire(Pg1Bls27,Pg2Bls27);
                        Print('Pairing e(aP,bQ): ',eBls27_2.toHexString,1);
                        if eBls27_1=eBls27_2 then Print('','Bilinearity Verified e(P,Q)^(a*b)= e(aP,bQ)' ,2)
                        else Print('','Bilinearity not Verified, somthing wrong' ,2);
                        Print('','-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'-------------------------------------------------------------------------------------------------------',0);
                        end
else if SelectedCurve=3 then begin
                        Print('Kss16 Curves Bilinearity Test :('+TreeView1.Selected.Text+' Curve.)' ,'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'----------------------------------------------------------------------------------------------------------------------------------',0);
                        Pg1Kss16:=KSS16CurvePairing1.HashToG1Point('kamel');
                        Pg2Kss16:=KSS16CurvePairing1.HashToG2Point('kamel');
                        Print('First Point P (from E(Fp)):',Pg1Kss16.toHexString,0);
                        Print('Second Point Q (from E(Fp2)):',Pg2Kss16.ToHexString,0);
                        tim:=now;
                        eKss16_1:=KSS16CurvePairing1.Paire(Pg1Kss16,Pg2Kss16);
                        tim:=now-tim;
                        DecodeTime(tim,h,m,s,ms);
                        time:='(in '+inttostr(m)+':'+inttostr(s)+':'+inttostr(ms)+'ms)';
                        Print('Pairing e(P,Q): '+time,eKss16_1.toHexString,0);
                        a.SetToRandom(KSS16CurvePairing1.P);
                        Print('First Scalar Exponent a (From Fp):',a.toHexString,0);
                        b.SetToRandom(KSS16CurvePairing1.P);
                        Print('Second Scalar Exponent b (From Fp):',b.toHexString,0);
                        ab:=a*b;
                        Print('Product a*b:',ab.toHexString,0);
                        eKss16_1:=eKss16_1.Pow(ab);
                        Print('Pairing''s exponent e(P,Q)^(a*b):',eKss16_1.toHexString,1);
                        Pg1Kss16:=a*Pg1Kss16;
                        Pg2Kss16:=b*Pg2Kss16;
                        Print('Point a*P :',Pg1Kss16.toHexString,7);
                        Print('Point b*Q :',Pg2Kss16.toHexString,8);
                        eKss16_2:=KSS16CurvePairing1.Paire(Pg1Kss16,Pg2Kss16);
                        Print('Pairing e(aP,bQ): ',eKss16_2.toHexString,1);
                        if eKss16_1=eKss16_2 then Print('','Bilinearity Verified e(P,Q)^(a*b)= e(aP,bQ)' ,2)
                        else Print('','Bilinearity not Verified, somthing wrong' ,2);
                        Print('','-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'-------------------------------------------------------------------------------------------------------',0);
                        end
else if SelectedCurve=4 then begin
                        Print('Kss18 Curves Bilinearity Test :('+TreeView1.Selected.Text+' Curve.)' ,'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'----------------------------------------------------------------------------------------------------------------------------------',0);
                        Pg1Kss18:=Kss18CurvePairing1.HashToG1Point('kamel');
                        Pg2Kss18:=Kss18CurvePairing1.HashToG2Point('kamel');
                        Print('First Point P (from E(Fp)):',Pg1Kss18.toHexString,0);
                        Print('Second Point Q (from E(Fp2)):',Pg2Kss18.ToHexString,0);
                        tim:=now;
                        eKss18_1:=Kss18CurvePairing1.Paire(Pg1Kss18,Pg2Kss18);
                        tim:=now-tim;
                        DecodeTime(tim,h,m,s,ms);
                        time:='(in '+inttostr(m)+':'+inttostr(s)+':'+inttostr(ms)+'ms)';
                        Print('Pairing e(P,Q): '+time,eKss18_1.toHexString,0);
                        a.SetToRandom(Kss18CurvePairing1.P);
                        Print('First Scalar Exponent a (From Fp):',a.toHexString,0);
                        b.SetToRandom(Kss18CurvePairing1.P);
                        Print('Second Scalar Exponent b (From Fp):',b.toHexString,0);
                        ab:=a*b;
                        Print('Product a*b:',ab.toHexString,0);
                        eKss18_1:=eKss18_1.Pow(ab);
                        Print('Pairing''s exponent e(P,Q)^(a*b):',eKss18_1.toHexString,1);
                        Pg1Kss18:=a*Pg1Kss18;
                        Pg2Kss18:=b*Pg2Kss18;
                        Print('Point a*P :',Pg1Kss18.toHexString,7);
                        Print('Point b*Q :',Pg2Kss18.toHexString,8);
                        eKss18_2:=Kss18CurvePairing1.Paire(Pg1Kss18,Pg2Kss18);
                        Print('Pairing e(aP,bQ): ',eKss18_2.toHexString,1);
                        if eKss18_1=eKss18_2 then Print('','Bilinearity Verified e(P,Q)^(a*b)= e(aP,bQ)' ,2)
                        else Print('','Bilinearity not Verified, somthing wrong' ,2);
                        Print('','-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'-------------------------------------------------------------------------------------------------------',0);
                        end
else if SelectedCurve=5 then begin
                        Print('Kss36 Curves Bilinearity Test :('+TreeView1.Selected.Text+' Curve.)' ,'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'----------------------------------------------------------------------------------------------------------------------------------',0);
                        Pg1Kss36:=KSS36Curve1.HashToG1Point('kamel');
                        Pg2Kss36:=KSS36Curve1.HashToG2Point('kamel');
                        Print('First Point P (from E(Fp)):',Pg1Kss36.toHexString,0);
                        Print('Second Point Q (from E(Fp2)):',Pg2Kss36.ToHexString,0);
                        tim:=now;
                        eKss36_1:=KSS36Curve1.Paire(Pg1Kss36,Pg2Kss36);
                        tim:=now-tim;
                        DecodeTime(tim,h,m,s,ms);
                        time:='(in '+inttostr(m)+':'+inttostr(s)+':'+inttostr(ms)+'ms)';
                        Print('Pairing e(P,Q): '+time,eKss36_1.toHexString,0);
                        a.SetToRandom(KSS36Curve1.P);
                        Print('First Scalar Exponent a (From Fp):',a.toHexString,0);
                        b.SetToRandom(KSS36Curve1.P);
                        Print('Second Scalar Exponent b (From Fp):',b.toHexString,0);
                        ab:=a*b;
                        Print('Product a*b:',ab.toHexString,0);
                        eKss36_1:=eKss36_1.Pow(ab);
                        Print('Pairing''s exponent e(P,Q)^(a*b):',eKss36_1.toHexString,1);
                        Pg1Kss36:=a*Pg1Kss36;
                        Pg2Kss36:=b*Pg2Kss36;
                        Print('Point a*P :',Pg1Kss36.toHexString,7);
                        Print('Point b*Q :',Pg2Kss36.toHexString,8);
                        eKss36_2:=KSS36Curve1.Paire(Pg1Kss36,Pg2Kss36);
                        Print('Pairing e(aP,bQ): ',eKss36_2.toHexString,1);
                        if eKss36_1=eKss36_2 then Print('','Bilinearity Verified e(P,Q)^(a*b)= e(aP,bQ)' ,2)
                        else Print('','Bilinearity not Verified, somthing wrong' ,2);
                        Print('','-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'-------------------------------------------------------------------------------------------------------',0);
                        end
                        else if SelectedCurve=7 then begin
                        case ComboBox1.ItemIndex of
                            0:MNTCurve1.CoordinatesSystem:=csAffine;
                            1:MNTCurve1.CoordinatesSystem:=csProjective;
                        end;
                        case ComboBox3.ItemIndex of
                        0:MNTCurve1.LoopMode:=lpmBinary;
                        1:MNTCurve1.LoopMode:=lpmNaf;
                        end;
                        case ComboBox4.ItemIndex of
                            0:MNTCurve1.PairingAlgorithm:=mnTate;
                            1:MNTCurve1.PairingAlgorithm:=mnAte;
                        end;
                        Print('Mnt Curves Bilinearity Test :('+TreeView1.Selected.Text+' Curve.)' ,'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'----------------------------------------------------------------------------------------------------------------------------------',0);
                        Pg1Mnt:=MntCurve1.HashToG1Point('kamel');
                        Pg2Mnt:=MntCurve1.HashToG2Point('kamel');
                        Print('First Point P (from E(Fp)):',Pg1Mnt.toHexString,0);
                        Print('Second Point Q (from E(Fp2)):',Pg2Mnt.ToHexString,0);
                        tim:=now;
                        eMnt_1:=MntCurve1.Paire(Pg1Mnt,Pg2Mnt);
                        tim:=now-tim;
                        DecodeTime(tim,h,m,s,ms);
                        time:='(in '+inttostr(m)+':'+inttostr(s)+':'+inttostr(ms)+'ms)';
                        Print('Pairing e(P,Q): '+time,eMnt_1.toHexString,0);
                        a.SetToRandom(MntCurve1.P);
                        Print('First Scalar Exponent a (From Fp):',a.toHexString,0);
                        b.SetToRandom(MntCurve1.P);
                        Print('Second Scalar Exponent b (From Fp):',b.toHexString,0);
                        ab:=a*b;
                        Print('Product a*b:',ab.toHexString,0);
                        eMnt_1:=eMnt_1.Pow(ab);
                        Print('Pairing''s exponent e(P,Q)^(a*b):',eMnt_1.toHexString,1);
                        Pg1Mnt:=a*Pg1Mnt;
                        Pg2Mnt:=b*Pg2Mnt;
                        Print('Point a*P :',Pg1Mnt.toHexString,7);
                        Print('Point b*Q :',Pg2Mnt.toHexString,8);
                        eMnt_2:=MntCurve1.Paire(Pg1Mnt,Pg2Mnt);
                        Print('Pairing e(aP,bQ): ',eMnt_2.toHexString,1);
                        if eMnt_1=eMnt_2 then Print('','Bilinearity Verified e(P,Q)^(a*b)= e(aP,bQ)' ,2)
                        else Print('','Bilinearity not Verified, somthing wrong' ,2);
                        Print('','-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'-------------------------------------------------------------------------------------------------------',0);
                        end
else if SelectedCurve=8 then begin
                        Print('Mnt Curves Bilinearity Test :('+TreeView1.Selected.Text+' Curve.)' ,'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'----------------------------------------------------------------------------------------------------------------------------------',0);
                        Pg1Mnt:=MntCurve1.HashToG1Point('kamel');
                        Pg2Mnt:=MntCurve1.HashToG2Point('kamel');
                        Print('First Point P (from E(Fp)):',Pg1Mnt.toHexString,0);
                        Print('Second Point Q (from E(Fp2)):',Pg2Mnt.ToHexString,0);
                        tim:=now;
                        eMnt_1:=MntCurve1.Paire(Pg1Mnt,Pg2Mnt);
                        tim:=now-tim;
                        DecodeTime(tim,h,m,s,ms);
                        time:='(in '+inttostr(m)+':'+inttostr(s)+':'+inttostr(ms)+'ms)';
                        Print('Pairing e(P,Q): '+time,eMnt_1.toHexString,0);
                        a.SetToRandom(MntCurve1.P);
                        Print('First Scalar Exponent a (From Fp):',a.toHexString,0);
                        b.SetToRandom(MntCurve1.P);
                        Print('Second Scalar Exponent b (From Fp):',b.toHexString,0);
                        ab:=a*b;
                        Print('Product a*b:',ab.toHexString,0);
                        eMnt_1:=eMnt_1.Pow(ab);
                        Print('Pairing''s exponent e(P,Q)^(a*b):',eMnt_1.toHexString,1);
                        Pg1Mnt:=a*Pg1Mnt;
                        Pg2Mnt:=b*Pg2Mnt;
                        Print('Point a*P :',Pg1Mnt.toHexString,7);
                        Print('Point b*Q :',Pg2Mnt.toHexString,8);
                        eMnt_2:=MntCurve1.Paire(Pg1Mnt,Pg2Mnt);
                        Print('Pairing e(aP,bQ): ',eMnt_2.toHexString,1);
                        if eMnt_1=eMnt_2 then Print('','Bilinearity Verified e(P,Q)^(a*b)= e(aP,bQ)' ,2)
                        else Print('','Bilinearity not Verified, somthing wrong' ,2);
                        Print('','-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'+'-------------------------------------------------------------------------------------------------------',0);
                        end;

end;

procedure TForm2.FormShow(Sender: TObject);
begin
TreeView1.FullExpand;
TreeView1.Selected:=TreeView1.Items[0];
end;



procedure TForm2.Memo1Change(Sender: TObject);
begin
SendMessage(Memo1.handle, WM_VSCROLL, SB_BOTTOM, 0);
end;

procedure TForm2.TreeView11CustomDrawItem(Sender: TCustomTreeView;
  Node: TTreeNode; State: TCustomDrawState; var DefaultDraw: Boolean);
begin
  with Sender as TCustomTreeView do
  begin
    if node.Level=0 then begin
                         Canvas.Font.Color := clBlue;
                         Canvas.Font.Style:=Canvas.Font.Style+[fsbold];
                         end;
      if node.Text[1]=''  then Canvas.Font.Name:='Symbol';
      if Node.ImageIndex=2 then begin
                                Canvas.Font.Color:=clred;
                                Canvas.Font.Style:=Canvas.Font.Style+[fsbold];
                                end;


  end;
end;

procedure TForm2.TreeView1Click(Sender: TObject);
var node:TTreeNode;
begin
node:=TreeView1.Selected;
GroupBox1.Enabled:=true;
if (node.Level=2)or(node.Parent.Index=0)
  or(node.Parent.Index=3)or(node.Parent.Index=4)then
  begin
  if (node.Level=2) then begin
                         case Node.Parent.Parent.Index of
                         1:begin
                           case node.Parent.Index of
                            0:begin
                              case node.Index of
                              0:BLS12CurvePairing1.Parametres:=sc128_1;
                              1:BLS12CurvePairing1.Parametres:=sc128_3;
                              2:BLS12CurvePairing1.Parametres:=sc128_4;
                              3:BLS12CurvePairing1.Parametres:=sc128_2_raz;
                              end;
                              BLS12CurvePairing1.GenerateParamsTreeView(TreeView11);
                              BLS12Parms;
                              SelectedCurve:=0;
                              end;
                            1:begin
                              case node.Index of
                              0:BLS24CurvePairing1.Parametres:=scBLS24at256_0;
                              1:BLS24CurvePairing1.Parametres:=scBLS24at256_1;
                              2:BLS24CurvePairing1.Parametres:=scBLS24at256_2;
                              3:BLS24CurvePairing1.Parametres:=scBLS24Razat256;
                              4:BLS24CurvePairing1.Parametres:=scBLS24at192_1;
                              5:BLS24CurvePairing1.Parametres:=scBLS24at192_2;
                              6:BLS24CurvePairing1.Parametres:=scBLS24at192_3;
                              7:BLS24CurvePairing1.Parametres:=scBLS24at192_Raz;
                              8:BLS24CurvePairing1.Parametres:=scBLS24at320;
                              end;
                              BLS24CurvePairing1.GenerateParamsTreeView(TreeView11);
                              BLS24Parms;
                              SelectedCurve:=1;
                              end;
                            2:begin
                              BLS27Curve1.GenerateParamsTreeView(TreeView11);
                              BLS27Parms;
                              SelectedCurve:=2;
                              end;
                           end;
                           end;
                         2:begin
                           case node.Parent.Index of
                            0:begin
                              case node.Index of
                              0:KSS16CurvePairing1.Parametres:=sc128Kss16_1;
                              1:KSS16CurvePairing1.Parametres:=sc128Kss16_3;
                              2:KSS16CurvePairing1.Parametres:=sc128Kss16_2_raz;
                              3:KSS16CurvePairing1.Parametres:=sc192Kss16_1;
                              4:KSS16CurvePairing1.Parametres:=sc192Kss16_2;
                              end;
                              KSS16CurvePairing1.GenerateParamsTreeView(TreeView11);
                              KSS16Parms;
                              SelectedCurve:=3;
                              end;
                            1:begin
                              case node.index of
                              0:KSS18CurvePairing1.Parametres:=scKSS18at192_1;
                              1:KSS18CurvePairing1.Parametres:=scKSS18at192_3;
                              2:KSS18CurvePairing1.Parametres:=scKSS18at192_2;
                              3:KSS18CurvePairing1.Parametres:=scKSS18at128;
                              4:KSS18CurvePairing1.Parametres:=scKSS18at256;
                              end;
                              KSS18CurvePairing1.GenerateParamsTreeView(TreeView11);
                              KSS18params;
                              SelectedCurve:=4;
                              end;
                            2:begin
                              case node.Index of
                              0:KSS36Curve1.Parametres:=scKSS36at192_1;
                              1:KSS36Curve1.Parametres:=scKSS36at192_2;
                              2:KSS36Curve1.Parametres:=scKSS36at256_1;
                              3:KSS36Curve1.Parametres:=scKSS36at256_2;
                              end;
                              KSS36Curve1.GenerateParamsTreeView(TreeView11);
                              kss36params;
                              SelectedCurve:=5;
                            end;
                           end;
                           end;
                         end;
                         end
  else begin
       case Node.Parent.Index of
       0:begin
         case Node.Index of
         0:BNCurvePairing1.Parametres:=scBeuchat;
         1:BNCurvePairing1.Parametres:=scAranha;
         2:BNCurvePairing1.Parametres:=scScott;
         3:BNCurvePairing1.Parametres:=scBCMNPZ;
         4:BNCurvePairing1.Parametres:=scISO224;
         5:BNCurvePairing1.Parametres:=scISO256;
         6:BNCurvePairing1.Parametres:=scISO384;
         7:BNCurvePairing1.Parametres:=scISO512;
         8:BNCurvePairing1.Parametres:=scRazvan;
         end;
         BNCurvePairing1.GenerateParamsTreeView(TreeView11);
         BNParms;
         SelectedCurve:=6;
         end;
       3:begin
         case node.Index of
         0:MNTCurve1.Parametres:=scCurveA;
         1:MNTCurve1.Parametres:=scCurveB;
         2:MNTCurve1.Parametres:=scCurveC;
         end;
         MNTCurve1.GenerateParamsTreeView(TreeView11);
         Mntparams;
         SelectedCurve:=7;
         end;
       4:begin
         case node.Index of
         0:SuperSingularCurvePairing1.Parametres:=scSS40;
         1:SuperSingularCurvePairing1.Parametres:=scSS66;
         2:SuperSingularCurvePairing1.Parametres:=scSS90;
         end;
         SuperSingularCurvePairing1.GenerateParamsTreeView(TreeView11);
         SSParms;
         SelectedCurve:=8;
         end;
       end;
       end;
  end
  else begin
       SelectedCurve:=-1;
       GroupBox1.Enabled:=false;
       end;
end;

end.
