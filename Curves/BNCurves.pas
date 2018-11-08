unit BNCurves;

interface
uses  vcl.dialogs,Treeview1,Vcl.ComCtrls, System.SysUtils,LargeIntegers,Fp12Arithmetic,Fp6Arithmetic,Fp2Arithmetic,GeneralTypes, System.classes,Vcl.forms,ECCFp2,ECCFp;


Type
   G1BN=FpPoint;
   G2BN=Fp2Point;
   GTBN=Fp12Int;
   LargeInt=Lint;
   BnCurvesPairingAlgos=(bpOptAte,bpR_Ate,bpEta,bpTate);
   StandardBNCurves=(scBeuchat,scAranha,scScott,scBCMNPZ,scISO224,scISO256,scISO384,scISO512,scRazvan);

   BNCurvesParamsDefinition=record
                            Identifier:String;
                            u:String;  // the paramater of generation for the BN curve
                            Beta:Integer; // non-square elements of the irredictible polynomial on FP2
                            Sigma:string; // non-square non-cube elements of the irredictible polynomial on FP6
                            A,B:String; // parametres of the curve   y^2=x^3+Ax+B
                            TwistMode:TTwistModel;
                            GenratorX:String;
                            TwistGeneratorBasePointX:String;
                            TwistGeneratorSeed:Word;
                            SecurityLevel:Integer;
                            end;
   ListOfBNParams=array of BNCurvesParamsDefinition;

 TBNCurvePairing=class(TCurve)
                      {Definition of a Curve with respect the the Weistrass model   Y^2=X^3+A*X+B (mod P)}
              private
              _LoopMode:LoopPoweringMode;                     //  Loopmode :Negative/Binary Representation
              _Fp12PoweringMode:FpPoweringMode;             // Final Exponentiation Powering Mode : Beuchat Approach, Karabina Approach, Classical Square and Multiply
              _PairingAlgorithm:BnCurvesPairingAlgos;         //  Pairing Algorithme OptAte,R-Ate.....
              Loop:PLIntArrayForm;
              _A,_B:Lint;
              _u:Lint;
              params:StandardBNCurves;
              procedure SetA(Value:String);
              function getA:String;
              procedure SetB(Value:String);
              function getB:String;
              procedure Setu(Value:String);
              function getu:String;
              function getRtw:String;
              function getAtw:String;
              function getBtw:String;
              procedure SetSigma(Value:String);
              function getSigma:String;
              procedure SetBeta(Value:integer);
              function getBeta:integer;
              function getR:string;
              function getN:string;
              function getTr:string;
              function getLp:string;
              function getHtw:string;
              function getH:string;
              function getP:string;
              procedure setCoord(value:CoordinatesSystem);
              function getCoord:CoordinatesSystem;
              function GetTwm:TTwistModel;
              procedure RecomputeParametres;
              procedure SetPArams(const Value: StandardBNCurves);
              Procedure FinalPowerOptimalAteBN(f:Fp12Int;u:LInt;PoweringMode:FpPoweringMode;var Result:Fp12Int); // Compute f^((p^12-1)/r)
              Function OptAtePairing(Pt:FpPoint;Qt:Fp2Point):Fp12Int;  // Compute Optimal Ate Pairings
              Function R_AtePairing(Pt:FpPoint;Qt:Fp2Point):Fp12Int; // Compute R-Ate Pairings
              Function EtaPairing(Pt:FpPoint;Qt:Fp2Point):Fp12Int;   // Compute Eta Pairings
              Function TatePairing(Pt:FpPoint;Qt:Fp2Point):Fp12Int;  // Compute Tate Pairings
              procedure SetLoopMode(Value:LoopPoweringMode);

            public
              Identifier:String;    /// A symbolic identifier of the curve
              property Beta:integer read getBeta write SetBeta;
              property Sigma:string read Getsigma write SetSigma;
              property u :string read Getu write Setu;
              property A:string read getA write setA;
              property B:string read getB write setB;
              property Atw :string read GetAtw;
              property Btw :string read GetBtw;
              property P:string read GetP;
              property H:string read GetH;
              property N:string read getN;
              property R:String read getR;
              property Rtw:string read getRtw;
              property Htw:string read getHtw;
              property Tr:string read getTr;
              property Lp:string read getLp;
              property TwistMode:TTwistModel read GetTwm;
              procedure SetStandardCurveParametres(inParametres:StandardBNCurves);
              procedure SetCustomCurveParametres(inParametres:BNCurvesParamsDefinition);

              constructor Create(AOwner : TComponent); override;
              destructor Destroy;
              function Paire(P:FpPoint;Q:Fp2Point):Fp12Int;
              function GetRandomG1Point:FpPoint;
              function GetRandomG2Point:Fp2Point;
              function GetDefautG1Generator:FpPoint;
              function GetDefautG2Generator:Fp2Point;
              function HashToG1Point(id:string):FpPoint;
              function HashToG2Point(id:String):Fp2Point;
              procedure GenerateParamsTreeView(Tree:TTreeView);override;
              published
              property CoordinatesSystem:CoordinatesSystem read getCoord write SetCoord;
              property Fp12PoweringMode:FpPoweringMode read _Fp12PoweringMode write _Fp12PoweringMode;
              property PairingAlgorithm:BnCurvesPairingAlgos read _PairingAlgorithm write _PairingAlgorithm ;
              property LoopMode:LoopPoweringMode read _LoopMode write SetLoopMode;
              property Parametres:StandardBNCurves read params write Setparams;
            end;




const
      BN_Beuchat254_Params:BNCurvesParamsDefinition=(Identifier:'Beuchat254';u:'0x3fc0100000000000';Beta:-5;sigma:'0+u*1';A:'0';B:'5';TwistMode:twDType;GenratorX:'-1';TwistGeneratorBasePointX:'1+u*0');
      BN_Aranha254_Params:BNCurvesParamsDefinition=(Identifier:'Aranha254';u:'-0x4080000000000001';Beta:-1;sigma:'1+u*1';A:'0';B:'2';TwistMode:twDType;GenratorX:'-1';TwistGeneratorBasePointX:'0+u*-1');
      BN_Scott254_Params:BNCurvesParamsDefinition=(Identifier:'Scott254';u:'-0x4000806000004081';Beta:-1;sigma:'1+u*1';A:'0';B:'2';TwistMode:twDType;GenratorX:'-1';TwistGeneratorBasePointX:'0+u*-1');
      BN_BCMNPZ254_Params:BNCurvesParamsDefinition=(Identifier:'BCMNPZ254';u:'-0x4000020100608205';Beta:-1;sigma:'1+u*1';A:'0';B:'2';TwistMode:twDType;GenratorX:'-1';TwistGeneratorBasePointX:'0+u*-1');
      BN_ISO224_Params:BNCurvesParamsDefinition=(Identifier:'ISO224';u:'-27162335522062337';Beta:-1;sigma:'1+u*1';A:'0';B:'3';TwistMode:twMType;GenratorX:'1';TwistGeneratorBasePointX:'2+u*1');
      BN_ISO256_Params:BNCurvesParamsDefinition=(Identifier:'ISO256';u:'-6990713071341666305';Beta:-1;sigma:'1+u*1';A:'0';B:'3';TwistMode:twMType;GenratorX:'1';TwistGeneratorBasePointX:'1+u*0');
      BN_ISO384_Params:BNCurvesParamsDefinition=(Identifier:'ISO320';u:'-29710560942993241785654870017';Beta:-1;sigma:'1+u*1';A:'0';B:'3';TwistMode:twMType;GenratorX:'1';TwistGeneratorBasePointX:'1+u*0');
      BN_ISO512_Params:BNCurvesParamsDefinition=(Identifier:'ISO512';u:'-128935115601040359985952327046386418689';Beta:-1;sigma:'1+u*1';A:'0';B:'3';TwistMode:twMType;GenratorX:'1';TwistGeneratorBasePointX:'1+u*0');
      BN_Razvan_update:BNCurvesParamsDefinition=(Identifier:'Razvan128Level';u:'0x4001FFFFFFFFFFFFFFFFFFFFFBFFF';Beta:-1;sigma:'1+u*1';A:'0';B:'-4';TwistMode:twMType;GenratorX:'2';TwistGeneratorBasePointX:'1+u*2');

procedure Register;

implementation


    {Generate parameters of a BN curves with precomputed constants from a compact definition }

Procedure FastLoadParams(param:StandardBNCurves;var Result:PtrCurveParams);
var rho,tmp,e:Lint;
    i:integer;
begin
if Result=nil then begin
                   new(Result);
                   new(Result.LoopBin);
                   new(Result.LoopNaf);
                   new(Result.LoopRateBin);
                   new(Result.LoopRateNaf);
                   new(Result.LoopEtaNaf);
                   new(Result.LoopEtaBin);
                   new(Result.LoopTateNaf);
                   new(Result.LoopTateBin);
                   end;
with Result^ do begin
case param of
scAranha:begin
               SecurityLevel:='Aranha254';
               Family:=cfBN;
               u:='-0x4080000000000001';
               P:='0x2523648240000001BA344D80000000086121000000000013A700000000000013';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0x2523648240000001BA344D8000000007FF9F800000000010A10000000000000D';
               Tr:='0x61818000000000030600000000000007';
               Lp:= '0x18300000000000004';
               n:='0x2523648240000001BA344D8000000007FF9F800000000010A10000000000000D';
               rho:='0x9366C48000000005B696800000000013A700000000000017';
               htw:='0x2523648240000001BA344D8000000008C2A2800000000016AD00000000000019';
               rtw:=r;
               A:=0;
               B:=2;
               H:=1;
               FieldParam:=InitFieldParams(-1,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam.Sigma.SetFormString('1+u*1');
               Atw.SetFormString('0');
               TowerParam^.FieldParam:=FieldParam;
               Atw.SetFieldParams(TowerParam^.FieldParam);
               Atw.SetFormString('0');
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               Btw.SetFieldParams(TowerParam^.FieldParam);
               TwistMode:=twDType;
               Btw:=B*TowerParam.Sigma.Inverse;
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopRateBin^:=(tr-1).ToIntArray;
               LoopRateNaf^:=(tr-1).ToNafArray;
               LoopEtaBin^:=rho.ToIntArray;
               LoopEtaNaf^:=rho.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointX:='-1';
               BasePointY:='0x2523648240000001BA344D80000000086121000000000013A700000000000012';
               TwistBasePointX.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointy.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointX.SetFormString('0+u*-1');
               TwistBasePointY.SetFormString('1');
        end;
scBeuchat:begin
               SecurityLevel:='Beuchat254';
               Family:=cfBN;
               u:='0x3FC0100000000000';
               P:='0x2370FB049D410FBE4E761A9886E502417D023F40180000017E80600000000001';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0x2370FB049D410FBE4E761A9886E502411DC1AF70120000017E80600000000001';
               Tr:='0x5F408FD0060000000000000000000001';
               Lp:= '0x17E80600000000002';
               n:='0x2370FB049D410FBE4E761A9886E502411DC1AF70120000017E80600000000001';
               rho:='0x8E521A9886E502411DC1AF70120000017E80600000000001';
               htw:='0x2370FB049D410FBE4E761A9886E50241DC42CF101E0000017E80600000000001';
               rtw:=r;
               A:=0;
               B:=5;
               H:=1;
               FieldParam:=InitFieldParams(-5,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam.Sigma.SetFormString('0+u*1');
               Atw.SetFormString('0');
               TowerParam^.FieldParam:=FieldParam;
               Atw.SetFieldParams(TowerParam^.FieldParam);
               Atw.SetFormString('0');
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               Btw.SetFieldParams(TowerParam^.FieldParam);
               TwistMode:=twDType;
               Btw:=B*TowerParam.Sigma.Inverse;
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopRateBin^:=(tr-1).ToIntArray;
               LoopRateNaf^:=(tr-1).ToNafArray;
               LoopEtaBin^:=rho.ToIntArray;
               LoopEtaNaf^:=rho.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointX:='-1';
               BasePointY:='0x2370FB049D410FBE4E761A9886E502417D023F40180000017E805FFFFFFFFFFF';
               TwistBasePointX.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointy.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointX.SetFormString('-1');
               TwistBasePointY.SetFormString('0xE86A09616F21FF57D258D5A8383F5ECB1176627C7036498A09DFF44CC78DEA9+u*0xB6097CB80DCE231E9E5E90DED770356AE84756E7D41259B244DC8A674D16996');
          end;
scScott:begin
               SecurityLevel:='Scott254';
               Family:=cfBN;
               u:='-0x4000806000004081';
               P:='0x240120DB6517014EFA0BAB3696F8D5F06E8A555614F464BABE9DBBFEEEB4A713';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0x240120DB6517014EFA0BAB3696F8D5F00E88D43492B2CB363A75777E8D30210D';
               Tr:='0x60018121824199848428448061848607';
               Lp:='0x18003024000018304';
               n:='0x240120DB6517014EFA0BAB3696F8D5F00E88D43492B2CB363A75777E8D30210D';
               rho:='0x9003628ECA2A09940C9D88B06C3A2EC9E84970AC7B2A2717';
               htw:='0x240120DB6517014EFA0BAB3696F8D5F0CE8BD6779735FE3F42C6007F50392D19';
               rtw:=r;
               A:=0;
               B:=2;
               H:=1;
               FieldParam:=InitFieldParams(-1,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam.Sigma.SetFormString('1+u*1');
               Atw.SetFormString('0');
               TowerParam^.FieldParam:=FieldParam;
               Atw.SetFieldParams(TowerParam^.FieldParam);
               Atw.SetFormString('0');
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               Btw.SetFieldParams(TowerParam^.FieldParam);
               TwistMode:=twDType;
               Btw:=B*TowerParam.Sigma.Inverse;
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopRateBin^:=(tr-1).ToIntArray;
               LoopRateNaf^:=(tr-1).ToNafArray;
               LoopEtaBin^:=rho.ToIntArray;
               LoopEtaNaf^:=rho.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointX:='-1';
               BasePointY:='0x240120DB6517014EFA0BAB3696F8D5F06E8A555614F464BABE9DBBFEEEB4A712';
               TwistBasePointX.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointy.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointX.SetFormString('0+u*-1');
               TwistBasePointY.SetFormString('1');
          end;
scBCMNPZ:begin
               SecurityLevel:='BCMNPZ254';
               Family:=cfBN;
               u:='-0x4000020100608205';
               P:='0x24000482410F5AADB74E200F3B89D00081CF93E428F0D651E8B2DC2BB460A48B';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0x24000482410F5AADB74E200F3B89D00021CF8DE127B73833D7FB71A511AA2BF5';
               Tr:='0x6000060301399E1E10B76A86A2B67897';
               Lp:='0x180000C0602430C1C';
               n:='0x240120DB6517014EFA0BAB3696F8D5F00E88D43492B2CB363A75777E8D30210D';
               rho:='0x90000D86C2F7D9E58CEAC8F763D9687528C67610746ACBEF';
               htw:='0x24000482410F5AADB74E200F3B89D000E1CF99E72A2A746FF96A46B257171D21';
               rtw:=r;
               A:=0;
               B:=2;
               H:=1;
               FieldParam:=InitFieldParams(-1,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam.Sigma.SetFormString('1+u*1');
               Atw.SetFormString('0');
               TowerParam^.FieldParam:=FieldParam;
               Atw.SetFieldParams(TowerParam^.FieldParam);
               Atw.SetFormString('0');
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               Btw.SetFieldParams(TowerParam^.FieldParam);
               TwistMode:=twDType;
               Btw:=B*TowerParam.Sigma.Inverse;
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopRateBin^:=(tr-1).ToIntArray;
               LoopRateNaf^:=(tr-1).ToNafArray;
               LoopEtaBin^:=rho.ToIntArray;
               LoopEtaNaf^:=rho.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointX:='-1';
               BasePointY:='0x24000482410F5AADB74E200F3B89D00081CF93E428F0D651E8B2DC2BB460A48A';
               TwistBasePointX.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointy.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointX.SetFormString('0+u*-1');
               TwistBasePointY.SetFormString('1');
          end;
scISO224: begin
               SecurityLevel:='ISO224';
               Family:=cfBN;
               u:='-27162335522062337';
               P:='0xBA139F3E23F1CDD79D9108355EA6E937BF36CC87626F8404E4E00013';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0xBA139F3E23F1CDD79D9108355EA60EF63EEE242757DD7E042420000D';
               Tr:='0xDA418048A8600A920600C0C00007';
               Lp:='0x243000060600004';
               n:='0xBA139F3E23F1CDD79D9108355EA60EF63EEE242757DD7E042420000D';
               rho:='0x1EDA225767F39E5D16B3C88A462B81B5A04E4E00017';
               htw:='0xBA139F3E23F1CDD79D9108355EA7C3793F7F74E76D018A05A5A00019';
               rtw:=r;
               A:=0;
               B:=3;
               H:=1;
               FieldParam:=InitFieldParams(-1,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam.Sigma.SetFormString('1+u*1');
               Atw.SetFormString('0');
               TowerParam^.FieldParam:=FieldParam;
               Atw.SetFieldParams(TowerParam^.FieldParam);
               Atw.SetFormString('0');
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               Btw.SetFieldParams(TowerParam^.FieldParam);
               TwistMode:=twMType;
               Btw:=B*TowerParam.Sigma;
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopRateBin^:=(tr-1).ToIntArray;
               LoopRateNaf^:=(tr-1).ToNafArray;
               LoopEtaBin^:=rho.ToIntArray;
               LoopEtaNaf^:=rho.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointX:='-1';
               BasePointY:='2';
               TwistBasePointX.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointy.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointX.SetFormString('2+u*1');
               TwistBasePointY.SetFormString('0xF62FA490E40C12E12CB7D2EC6D8B46D74DF010F746FCBB57FB007C8+u*0xA7C962F8530A796C7B3EB5D8A487847BE1CFB1E4C9DB3C7C4CE6092B');
          end;
scISO256: begin
               SecurityLevel:='ISO256';
               Family:=cfBN;
               u:='-6990713071341666305';
               P:='0xBE15F189A88462F19AF66402A838D25702945A461021001D8F38270000000013';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0xBE15F189A88462F19AF66402A838D25625FC279FF81F8019030821000000000D';
               Tr:='0xDC9832A6180180048C30060000000007';
               Lp:='0x24618030000000004';
               n:='0xBE15F189A88462F19AF66402A838D25625FC279FF81F8019030821000000000D';
               rho:='0x1F596B40022BD9B78F16AF7BB6816801D8F38270000000017';
               htw:='0xBE15F189A88462F19AF66402A838D257DF2C8CEC282280221B682D0000000019';
               rtw:=r;
               A:=0;
               B:=3;
               H:=1;
               FieldParam:=InitFieldParams(-1,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam.Sigma.SetFormString('1+u*1');
               Atw.SetFormString('0');
               TowerParam^.FieldParam:=FieldParam;
               Atw.SetFieldParams(TowerParam^.FieldParam);
               Atw.SetFormString('0');
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               Btw.SetFieldParams(TowerParam^.FieldParam);
               TwistMode:=twMType;
               Btw:=B*TowerParam.Sigma;
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopRateBin^:=(tr-1).ToIntArray;
               LoopRateNaf^:=(tr-1).ToNafArray;
               LoopEtaBin^:=rho.ToIntArray;
               LoopEtaNaf^:=rho.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointX:='1';
               BasePointY:='2';
               TwistBasePointX.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointy.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointX.SetFormString('1+u*0');
               TwistBasePointY.SetFormString('0xC0A9DD63B67535E3C33E93B1464D645799BBE954F2935802032EB902114D1DF+u*0x99F61806F64E68D6E65AA8516B0A4F8695C11E8622A55F9D2E9F644F9CC18A76');
          end;
scISO384:begin
               SecurityLevel:='ISO384';
               Family:=cfBN;
               u:='-29710560942993241785654870017';
               P:='0xB64000000F300000007D4C05B201ECC05B20038B99FCC81521045300C6DEC2A66241A41DE2010802DC0D802100270013';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0xB64000000F300000007D4C05B201ECC05B20038B99FCC81449045300BDDEC2A6622764196200FC02C40D801F8021000D';
               Tr:='0xD800000009000000001A400480000C001800000180060007';
               Lp:='0x2400000000C00000000030004';
               n:='0xB64000000F300000007D4C05B201ECC05B20038B99FCC81449045300BDDEC2A6622764196200FC02C40D801F8021000';
               rho:='0x1E60000001E60000000A9980CA8017100870000E22189C01D7600B4009C04801680270017';
               htw:='0xB64000000F300000007D4C05B201ECC05B20038B99FCC815F9045300CFDEC2A6625BE42262011402F40D8022802D0019';
               rtw:=r;
               A:=0;
               B:=3;
               H:=1;
               FieldParam:=InitFieldParams(-1,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam.Sigma.SetFormString('1+u*1');
               Atw.SetFormString('0');
               TowerParam^.FieldParam:=FieldParam;
               Atw.SetFieldParams(TowerParam^.FieldParam);
               Atw.SetFormString('0');
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               Btw.SetFieldParams(TowerParam^.FieldParam);
               TwistMode:=twMType;
               Btw:=B*TowerParam.Sigma;
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopRateBin^:=(tr-1).ToIntArray;
               LoopRateNaf^:=(tr-1).ToNafArray;
               LoopEtaBin^:=rho.ToIntArray;
               LoopEtaNaf^:=rho.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointX:='1';
               BasePointY:='2';
               TwistBasePointX.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointy.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointX.SetFormString('1+u*0');
               TwistBasePointY.SetFormString('0x2705D177D06F7ACFC4021860560EBCA2BD2F3805E714A3998804BD734CB469F36F0BC2B4EBD4E6A06D01D066843DA267+u*0x412E8B989DE18F90B47702E4AFD5B6D823925B79E4BEDD4888F61AA6E0C184CC151E5BFF1E82542195080EED736E18DE');
         end;
scISO512:Begin
               SecurityLevel:='ISO512';
               Family:=cfBN;
               u:='-128935115601040359985952327046386418689';
               P:='0xBDF69624FAAC52007C0B60001B481F5B6C8A0B08A8888A3CEB10012037F0E7A404CBDAD0813614421000000DABADD04210888409C00000000000241B08413813';
               new(p.InverseFordivision);
               p.Limit:=2*p.Data.i16[-2];
               p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
               p.id:='Y';
               r:='0xBDF69624FAAC52007C0B60001B481F5B6C8A0B08A8888A3CEB10012037F0E7A32845DACFEFB61441F800000DABADBE0D84887E08400000000000241B07E1080D';
               Tr:='0xDC8600009180000018000000000012348C000601800000000000000000603007';
               Lp:='0x246000000C00000000000000000001804';
               n:='0xBDF69624FAAC52007C0B60001B481F5B6C8A0B08A8888A3CEB10012037F0E7A32845DACFEFB61441F800000DABADBE0D84887E08400000000000241B07E1080D';
               rho:='0x1F558A401F02D8000A3B0000012003E129BDA28F4868006C168000000028FD0ED8ED85A09C00000000000000905A13817';
               htw:='0xBDF69624FAAC52007C0B60001B481F5B6C8A0B08A8888A3CEB10012037F0E7A4E151DAD112B614422800000DABADE2769C888A0B400000000000241B08A16819';
               rtw:=r;
               A:=0;
               B:=3;
               H:=1;
               FieldParam:=InitFieldParams(-1,p);
               New(FieldParam.MontgomeryData);
               InitMontgomeryStruct(P,FieldParam.MontgomeryData);
               New(TowerParam);
               TowerParam.Sigma.SetFormString('1+u*1');
               Atw.SetFormString('0');
               TowerParam^.FieldParam:=FieldParam;
               Atw.SetFieldParams(TowerParam^.FieldParam);
               Atw.SetFormString('0');
               TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
               Btw.SetFieldParams(TowerParam^.FieldParam);
               TwistMode:=twMType;
               Btw:=B*TowerParam.Sigma;
               LoopBin^:=Lp.ToIntArray;
               LoopNaf^:=Lp.ToNafArray;
               LoopRateBin^:=(tr-1).ToIntArray;
               LoopRateNaf^:=(tr-1).ToNafArray;
               LoopEtaBin^:=rho.ToIntArray;
               LoopEtaNaf^:=rho.ToNafArray;
               LoopTateBin^:=n.ToIntArray;
               LoopTateNaf^:=n.ToNafArray;
               BasePointX:='1';
               BasePointY:='2';
               TwistBasePointX.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointy.SetFieldParams(TowerParam^.FieldParam);
               TwistBasePointX.SetFormString('1+u*0');
               TwistBasePointY.SetFormString('0x3DE23849792F45229D4C43485F72DA62E7A5AD0502E0625CDDD5624AE89195CEEB1A9C4C81D56EE75703516A398052E655339EE38DCFEE56ADFB84EF22F75ED6'+'+u*0x44FED488F1E8298A4269626FCEF9032B59903F99FE76326518FDA3F7E3C2637437C05EAFBB5C78C0AF60BCEFF2CD78F10EDA75F169034FBF60D954D9F5B1B91');
         End;
end;
TowerParam^.FrobeniusP_Const[0].SetFieldParams(TowerParam^.FieldParam);
TowerParam^.FrobeniusP_Const[0]:=TowerParam^.Sigma.Pow((P-1)/6);
for i:=1 to 4 do begin
                 TowerParam^.FrobeniusP_Const[i].SetFieldParams(TowerParam^.FieldParam);
                 TowerParam^.FrobeniusP_Const[i]:=TowerParam^.FrobeniusP_Const[i-1]*TowerParam^.FrobeniusP_Const[0];
                 end;
for i:=0 to 4 do begin
                 TowerParam^.FrobeniusP2_Const[i].SetFieldParams(TowerParam^.FieldParam);
                 TowerParam^.FrobeniusP3_Const[i].SetFieldParams(TowerParam^.FieldParam);
                 TowerParam^.FrobeniusP2_Const[i]:=TowerParam^.FrobeniusP_Const[i]*TowerParam^.FrobeniusP_Const[i].Conjugate;
                 TowerParam^.FrobeniusP3_Const[i]:=TowerParam^.FrobeniusP_Const[i]*TowerParam^.FrobeniusP2_Const[i];
                 end;
for i:=0 to 2 do begin
                 e:=P.Pow(i+1);
                 FrobeniusMapConstX[i].SetFieldParams (TowerParam.FieldParam);
                 case TwistMode of
                                 twDType: begin
                                          FrobeniusMapConstX[i]:=TowerParam^.Sigma.Pow((e-1)/3);
                                          FrobeniusMapConstY[i]:=TowerParam^.Sigma.Pow((e-1)/2);
                                          end;
                                 twMType:begin
                                         FrobeniusMapConstX[i]:=(TowerParam^.Sigma.Pow((e-1)/3)).inverse;
                                         FrobeniusMapConstY[i]:=TowerParam^.Sigma.Pow((e-1)/2).inverse;
                                         end;
                            end;
                 end;
_Pow_LInt(P,power,4);
_Pow_LInt(P,tmp,2);
_Sub_LInt(power,tmp,power);
_Inc_LInt(power,1);
_Div_Mod_LInt(power,R,power,tmp);
end;
end;

procedure GenerateBNParametres(Params: BNCurvesParamsDefinition; var Result:PtrCurveParams);
var localP:LInt;
    u4,u3,u2,tmp,tmp1:LInt;
    i:integer;
    e:LInt;
    twtmp,twtmp1:Fp2Int;
    rho:Lint; // Size of The Miller Loop for Eta Pairings
begin
if Result=nil then begin
                   new(Result);
                   new(Result.LoopBin);
                   new(Result.LoopNaf);
                   new(Result.LoopRateBin);
                   new(Result.LoopRateNaf);
                   new(Result.LoopEtaNaf);
                   new(Result.LoopEtaBin);
                   new(Result.LoopTateNaf);
                   new(Result.LoopTateBin);
                   end;
with Result^ do begin
                 SecurityLevel:=Params.Identifier;
                Family:=cfBN;
                if (Pos('0x',Params.u)<>0)or(Pos('$',Params.u)<>0) then u:=_hex_To_LInt(Params.u)
                else u:=_str_To_LInt(Params.u);
                u2:=u.Sqr;
                u3:=u*u2;
                u4:=u2.Sqr;
                usqr:=u2;
                ucub:=u3;
                  LocalP:=36*u4+36*u3+24*u2+6*u+1;
                if not IsLIntPrime(LocalP) then raise Exception.Create('Paramétre t invalide pour la construction de la courbe BN...')
                else begin
                     ///   should test if beta and sigma are valide parametres
                     p:=Localp;
                     new(p.InverseFordivision);
                     p.Limit:=2*p.Data.i16[-2];
                     p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
                     p.id:='Y';
                     R:=LocalP-(6*u2);
                     Tr:=6*u2+1;
                     Lp:=Lint(6*u+2).Absolute;
                     n:=p+1-Tr;
                     end;
                A:=_Str_To_LInt(Params.A);
                B:=_Str_To_LInt(Params.B);
                Htw:=36*u4+36*u3+30*u2+6*u+1;
                Rho:=(36*u3+18*u2+6*u+1).Absolute; // Miller's Loop Size for the Eta pairing
                H:=1;
                FieldParam:=InitFieldParams(Params.Beta,p);
                New(FieldParam.MontgomeryData);
                InitMontgomeryStruct(P,FieldParam.MontgomeryData);
                e:=params.beta;
                if e.IsASqrMod(P,FieldParam.MontgomeryData) then Exception.Create(inttostr(Params.beta)+' is a quadratic Residue modulo P......');
                New(TowerParam);
                TowerParam.Sigma.SetFormString(Params.Sigma);
                Atw.SetFormString(Params.A);
                TowerParam^.FieldParam:=FieldParam;
                Atw.SetFieldParams(TowerParam^.FieldParam);
                Atw.SetFormString(Params.A);
                TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
                Btw.SetFieldParams(TowerParam^.FieldParam);
                TwistMode:=Params.TwistMode;
                if TwistMode=twMType then Btw:=B*TowerParam.Sigma
                else Btw:=B*TowerParam.Sigma.Inverse;
                Rtw:=R; ///  The true order is R*(2*p-R) , so we ca choos a subgourp of order R (the cofactor will be 2p-R) (https://eprint.iacr.org/2005/133.pdf)
                LoopBin^:=Lp.ToIntArray;
                LoopNaf^:=Lp.ToNafArray;
                LoopRateBin^:=(tr-1).ToIntArray;
                LoopRateNaf^:=(tr-1).ToNafArray;
                LoopEtaBin^:=rho.ToIntArray;
                LoopEtaNaf^:=rho.ToNafArray;
                LoopTateBin^:=n.ToIntArray;
                LoopTateNaf^:=n.ToNafArray;
                if (Pos('0x',Params.GenratorX)<>0)or(Pos('$',Params.GenratorX)<>0) then BasePointX:=_hex_To_LInt(Params.GenratorX)
                else BasePointX:=_Str_To_LInt(Params.GenratorX);
                TwistBasePointX.SetFieldParams(TowerParam^.FieldParam);
                TwistBasePointy.SetFieldParams(TowerParam^.FieldParam);
                TwistBasePointX.SetFormString(Params.TwistGeneratorBasePointX);
                ///   Constructing the generator for G1
                _Sqr_LInt(BasePointX,tmp);
                _Add_LInt(tmp,A,tmp);
                _Mul_LInt(tmp,BasePointX,tmp1);
                _Add_LInt(tmp1,B,tmp);
                _Mod_LInt(tmp,P,tmp);
                if tmp.IsASqrMod(P,FieldParam.MontgomeryData)
                                              then begin
                                                   ModSquareLInt(tmp,P, BasePointY,FieldParam.MontgomeryData,false);
                                                   BasePointY:=P-BasePointY;
                                                   end
                 else raise Exception.Create('Invalide parameter for the generator of the curve');
                ///   Constructing the generator for G2
                _Sqr_FP2(TwistBasePointX,twtmp);
                _Add_FP2(twtmp,Atw,twtmp);
                _Mul_FP2(twtmp,TwistBasePointX,twtmp1);
                _Add_FP2(twtmp1,Btw,twtmp);
                if twtmp.IsASquare then begin
                                        _Sqrt_FP2(twtmp,TwistBasePointy);
                                        _Neg_FP2(TwistBasePointy,twtmp1);
                                        if _Compare_FP2(TwistBasePointy,twtmp1)=1 then TwistBasePointy:=twtmp1;
                                        end
                else raise Exception.Create('Invalide parameter for the generator of the twisted curve');
                TowerParam^.FrobeniusP_Const[0].SetFieldParams(TowerParam^.FieldParam);
                TowerParam^.FrobeniusP_Const[0]:=TowerParam^.Sigma.Pow((P-1)/6);
                for i:=1 to 4 do begin
                                 TowerParam^.FrobeniusP_Const[i].SetFieldParams(TowerParam^.FieldParam);
                                 TowerParam^.FrobeniusP_Const[i]:=TowerParam^.FrobeniusP_Const[i-1]*TowerParam^.FrobeniusP_Const[0];
                                 end;
                for i:=0 to 4 do begin
                                 TowerParam^.FrobeniusP2_Const[i].SetFieldParams(TowerParam^.FieldParam);
                                 TowerParam^.FrobeniusP3_Const[i].SetFieldParams(TowerParam^.FieldParam);
                                 TowerParam^.FrobeniusP2_Const[i]:=TowerParam^.FrobeniusP_Const[i]*TowerParam^.FrobeniusP_Const[i].Conjugate;
                                 TowerParam^.FrobeniusP3_Const[i]:=TowerParam^.FrobeniusP_Const[i]*TowerParam^.FrobeniusP2_Const[i];
                                 end;
                for i:=0 to 2 do begin
                                 e:=P.Pow(i+1);
                                 FrobeniusMapConstX[i].SetFieldParams (TowerParam.FieldParam);
                                 case TwistMode of
                                 twDType: begin
                                          FrobeniusMapConstX[i]:=TowerParam^.Sigma.Pow((e-1)/3);
                                          FrobeniusMapConstY[i]:=TowerParam^.Sigma.Pow((e-1)/2);
                                          end;
                                 twMType:begin
                                         FrobeniusMapConstX[i]:=(TowerParam^.Sigma.Pow((e-1)/3)).inverse;
                                         FrobeniusMapConstY[i]:=TowerParam^.Sigma.Pow((e-1)/2).inverse;
                                         end;
                                 end;
                                 end;
                _Pow_LInt(P,power,4);
                _Pow_LInt(P,tmp,2);
                _Sub_LInt(power,tmp,power);
                _Inc_LInt(power,1);
                _Div_Mod_LInt(power,R,power,tmp);
                end;
end;


{ TBNCurve }

{*******************************************************************************}
constructor TBNCurvePairing.Create(AOwner : TComponent);
begin
inherited create(AOwner);
CurveParams:=nil;
SetStandardCurveParametres(scAranha);
_A:=Curveparams.A;
_B:=Curveparams.B;
_u:=Curveparams.u;
CoordinatesSystem:=csProjective;
Fp12PoweringMode:=pmKarbina;
SetLoopMode(lpmAuto);
PairingAlgorithm:=bpOptAte;
end;

{*******************************************************************************}
destructor TBNCurvePairing.Destroy;
begin
Dispose(CurveParams);
inherited destroy;
end;

{*******************************************************************************}
function TBNCurvePairing.getA: String;
begin
Result:=_A.ToHexString;
end;

{*******************************************************************************}
function TBNCurvePairing.getAtw: String;
begin
Result:=CurveParams^.Atw.toHexString;
end;

{*******************************************************************************}
function TBNCurvePairing.getB: String;
begin
Result:=_B.ToHexString;
end;

{*******************************************************************************}
function TBNCurvePairing.getBeta:integer;
begin
Result:=CurveParams^.FieldParam.Beta;
end;

{*******************************************************************************}
function TBNCurvePairing.getBtw: String;
begin
Result:=CurveParams^.Btw.toHexString;
end;

{*******************************************************************************}
function TBNCurvePairing.getCoord: CoordinatesSystem;
begin
Result:=CurveParams^.CoordSys;
end;

{*******************************************************************************}
function TBNCurvePairing.GetDefautG1Generator: FpPoint;
begin
Result.SetToDefaultGenerator;
end;

{*******************************************************************************}
function TBNCurvePairing.GetDefautG2Generator: Fp2Point;
begin
Result.SetToDefaultGenerator;
end;

{*******************************************************************************}
function TBNCurvePairing.getH: string;
begin
Result:=CurveParams^.H.toHexString;
end;

{*******************************************************************************}
function TBNCurvePairing.getHtw: string;
begin
Result:=CurveParams^.Htw.toHexString;
end;

{*******************************************************************************}
function TBNCurvePairing.getLp: string;
begin
Result:=CurveParams^.Lp.toHexString;
end;

{*******************************************************************************}
function TBNCurvePairing.getN: string;
begin
Result:=CurveParams^.N.toHexString;
end;

{*******************************************************************************}
function TBNCurvePairing.getP: string;
begin
Result:=CurveParams^.P.toHexString;
end;

{*******************************************************************************}
function TBNCurvePairing.getR: string;
begin
Result:=CurveParams^.R.toHexString;
end;

{*******************************************************************************}
function TBNCurvePairing.GetRandomG1Point: FpPoint;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsRandomTorsionPoint;
end;

{*******************************************************************************}
function TBNCurvePairing.GetRandomG2Point: Fp2Point;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsRandomTorsionPoint;
end;

{*******************************************************************************}
function TBNCurvePairing.getRtw: String;
begin
Result:=CurveParams^.Rtw.toHexString;
end;

{*******************************************************************************}
function TBNCurvePairing.getSigma: String;
begin
Result:=CurveParams^.TowerParam.Sigma.toHexString;
end;

{*******************************************************************************}
function TBNCurvePairing.getTr: string;
begin
Result:=CurveParams^.Tr.toHexString;
end;

{*******************************************************************************}
function TBNCurvePairing.GetTwm: TTwistModel;
begin
Result:=CurveParams.TwistMode;
end;

{*******************************************************************************}
function TBNCurvePairing.getu: String;
begin
Result:=_u.ToHexString;
end;

function TBNCurvePairing.HashToG1Point(id: string): FpPoint;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsTorsionFromString(id);
end;

function TBNCurvePairing.HashToG2Point(id: String): Fp2Point;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsTorsionFromString(id);
end;

{*******************************************************************************}
procedure TBNCurvePairing.RecomputeParametres;
var tmp:BNCurvesParamsDefinition;
begin
tmp.u:=_u.toHexString;
tmp.Beta:=beta;
tmp.Sigma:=Sigma;
tmp.A:=_A.toHexString;
tmp.B:=_B.toHexString;
tmp.Identifier:=Identifier;
GenerateBNParametres(tmp,CurveParams);
end;

{*******************************************************************************}
procedure TBNCurvePairing.SetStandardCurveParametres(inParametres:StandardBNCurves);
begin
case inParametres of
    scBeuchat:begin
              FastLoadParams(scBeuchat,CurveParams);
              _A:=_Str_To_LInt(BN_Beuchat254_Params.A);
              _B:=_Str_To_LInt(BN_Beuchat254_Params.B);
              _u:=_Hex_To_LInt( BN_Beuchat254_Params.u);
              end;
    scAranha:begin
             FastLoadParams(scAranha,CurveParams);
             _A:=_Str_To_LInt(BN_Aranha254_Params.A);
             _B:=_Str_To_LInt(BN_Aranha254_Params.B);
             _u:=_Hex_To_LInt(BN_Aranha254_Params.u);
             end;
    scScott:begin
             FastLoadParams(scScott,CurveParams);
             _A:=_Str_To_LInt(BN_Scott254_Params.A);
             _B:=_Str_To_LInt(BN_Scott254_Params.B);
             _u:=_Hex_To_LInt(BN_Scott254_Params.u);
             end;
    scBCMNPZ:begin
             FastLoadParams(scBCMNPZ,CurveParams);
             _A:=_Str_To_LInt(BN_BCMNPZ254_Params.A);
             _B:=_Str_To_LInt(BN_BCMNPZ254_Params.B);
             _u:=_Hex_To_LInt(BN_BCMNPZ254_Params.u);
             end;
    scISO224:begin
             FastLoadParams(scISO224,CurveParams);
             _A:=_Str_To_LInt(BN_ISO224_Params.A);
             _B:=_Str_To_LInt(BN_ISO224_Params.B);
             _u:=_Hex_To_LInt(BN_ISO224_Params.u);
             end;
    scISO256:begin
             FastLoadParams(scISO256,CurveParams);
             _A:=_Str_To_LInt(BN_ISO256_Params.A);
             _B:=_Str_To_LInt(BN_ISO256_Params.B);
             _u:=_Hex_To_LInt(BN_ISO256_Params.u);
             end;
    scISO384:begin
             FastLoadParams(scISO384,CurveParams);
             _A:=_Str_To_LInt(BN_ISO384_Params.A);
             _B:=_Str_To_LInt(BN_ISO384_Params.B);
             _u:=_Hex_To_LInt(BN_ISO384_Params.u);
             end;
    scISO512:begin
             FastLoadParams(scISO512,CurveParams);
             _A:=_Str_To_LInt(BN_ISO512_Params.A);
             _B:=_Str_To_LInt(BN_ISO512_Params.B);
             _u:=_Hex_To_LInt(BN_ISO512_Params.u);
             end;
    scRazvan:begin
             GenerateBNParametres(BN_Razvan_update,CurveParams);
             _A:=_Str_To_LInt(BN_Razvan_update.A);
             _B:=_Str_To_LInt(BN_Razvan_update.B);
             _u:=_Hex_To_LInt(BN_Razvan_update.u);
             end;

end;
end;

{*******************************************************************************}
procedure TBNCurvePairing.SetA(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBNCurvePairing.SetB(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBNCurvePairing.SetBeta(Value: Integer);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBNCurvePairing.setCoord(value: CoordinatesSystem);
begin
CurveParams^.CoordSys:=value;
end;

{*******************************************************************************}
procedure TBNCurvePairing.SetCustomCurveParametres(inParametres: BNCurvesParamsDefinition);
begin
GenerateBNParametres(inParametres,CurveParams);
_A:=_Str_To_LInt(inParametres.A);
_B:=_Str_To_LInt(inParametres.B);
_u:=_Hex_To_LInt(inParametres.u);
end;

{*******************************************************************************}
procedure TBNCurvePairing.SetPArams(const Value: StandardBNCurves);
var pr:TComponent;
    i:integer;
begin
SetStandardCurveParametres(Value);
pr:=Self.Owner;
for i:=0 to pr.ComponentCount-1 do begin
                                   if (pr.Components[i] is TTreeView1)and (TTreeView1(pr.Components[i]).Curve = self)
                                   then GenerateParamsTreeView(TTreeView(pr.Components[i]));
                                   end;
params:=value;
end;

{*******************************************************************************}
procedure TBNCurvePairing.SetSigma(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBNCurvePairing.Setu(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBNCurvePairing.GenerateParamsTreeView(Tree:TTreeView);
var s:string;
    r:real;
    G:FpPoint;Gtw:Fp2Point;
    tmp:TTreeNode;
begin
Tree.Items.Clear;
tmp:=Tree.Items.Add(nil,'Curve');
Tree.Items.AddChild(tmp,'Family : BN ('+curveparams.SecurityLevel+')');
s:='y^2=x^3+'+CurveParams.B.ToDecimalString;
Tree.Items.AddChild(tmp,'Equation : '+s);
tmp:=Tree.Items.Add(nil,'Initial Parameter x0');
Tree.Items.AddChild(tmp,'x0= '+curveparams.u.ToHexString );
Tree.Items.AddChild(tmp,'Hamming Weight (Bin) ='+IntToStr(HammingWeight(CurveParams.u,False)));
Tree.Items.AddChild(tmp,'Hamming Weight (Neg) ='+IntToStr(HammingWeight(CurveParams.u,True)));
tmp:=Tree.Items.Add(nil,'Field Fp');
Tree.Items.AddChild(tmp,'Prime P = '+CurveParams.P.ToHexString );
Tree.Items.AddChild(tmp,'Size of Fp : '+inttostr(CurveParams.P.BitLength)+' bit');
tmp:=Tree.Items.Add(nil,'Field E(Fp)');
Tree.Items.AddChild(tmp,'Order (#E(Fp)) N = '+CurveParams.N.ToHexString );
Tree.Items.AddChild(tmp,'Frobenius Trace Tr : '+CurveParams.Tr.ToHexString);
tmp:=Tree.Items.Add(nil,'Base Fields G1 : E(Fp)[R]');
Tree.Items.AddChild(tmp,'Order of G1 R = '+CurveParams.R.ToHexString);
Tree.Items.AddChild(tmp,'Size of G1 : '+inttostr(CurveParams.R.BitLength)+' bit');
Tree.Items.AddChild(tmp,'Cofactor of G1 H = '+CurveParams.H.ToHexString);
tmp:=Tree.Items.Add(nil,'Twist');
Tree.Items.AddChild(tmp,'Equation of the Twist : '+'y^2=x^3+'+Curveparams.Btw.toHexString);
Tree.Items.AddChild(tmp,'Degree of the Twist : Sextic');
if Curveparams.TwistMode=twMType then S:='M-Type' else S:='D-Type';
Tree.Items.AddChild(tmp,'Type of the Twist : '+S);
tmp:=Tree.Items.Add(nil,'Field G2 (Twist) :E''(Fp2)[R]');
Tree.Items.AddChild(tmp,'Order of G2 R = '+CurveParams.Rtw.ToHexString);
Tree.Items.AddChild(tmp,'Size of G2 : '+inttostr(CurveParams.Rtw.BitLength)+' bit');
Tree.Items.AddChild(tmp,'Cofactor of G2 Htw = '+CurveParams.Htw.ToHexString);
tmp:=Tree.Items.Add(nil,'Targted Field (GT):');
Tree.Items.AddChild(tmp,'Extension Field : Fp12');
Tree.Items.AddChild(tmp,'Size of GT : '+inttostr(12*CurveParams.P.BitLength)+' bit');
tmp:=Tree.Items.Add(nil,'Security Level');
Tree.Items.AddChild(tmp,'Security in G1 : '+inttostr(CurveParams.R.BitLength div 2)+' bit').ImageIndex:=2;
Tree.Items.AddChild(tmp,'Security in G2 : '+inttostr(CurveParams.R.BitLength div 2)+' bit').ImageIndex:=2;
r:=1.92* exp((1/3)*ln(12*CurveParams.P.BitLength))*sqr(exp((1/3)*(ln(ln(12*CurveParams.P.BitLength)/ln(2)))));
Tree.Items.AddChild(tmp,'Security in GT (GNFS) : '+inttostr(Round(r))+' bit').ImageIndex:=2;
tmp:=Tree.Items.Add(nil,'Implemented Pairings');
Tree.Items.AddChild(tmp,'Tate');
Tree.Items.AddChild(tmp,'Otp-Ate');
Tree.Items.AddChild(tmp,'R-Ate');
tmp:=Tree.Items.Add(nil,'Tower Construction');
if Curveparams.FieldParam.Beta<0 then Tree.Items.AddChild(tmp,'Fp2<u>=ExstensionField<u,|u^2-'+inttostr(Abs(Curveparams.FieldParam.Beta))+'>')
else Tree.Items.AddChild(tmp,'Fp2<u>=ExstensionField<u,|u^2+'+inttostr(Curveparams.FieldParam.Beta)+'>');
Tree.Items.AddChild(tmp,'Fp6<v>=ExstensionField<v,|v^3-('+Curveparams.TowerParam.Sigma.toHexString+')>');
Tree.Items.AddChild(tmp,'Fp12<w>=ExstensionField<w,|w^2-v>');
tmp:=Tree.Items.Add(nil,' = '+floattostrf(CurveParams.P.BitLength/CurveParams.R.BitLength,ffGeneral, 2, 4 ));
tmp:=Tree.Items.Add(nil,'Default G1 Generator');
G.SetCurveParams (Self.CurveParams);
G.SetToDefaultGenerator;
Gtw.SetCurveParams (Self.CurveParams);
Gtw.SetToDefaultGenerator;
Tree.Items.AddChild(tmp,G.toHexString);
tmp:=Tree.Items.Add(nil,'Default G2 Generator');
Tree.Items.AddChild(tmp,Gtw.toHexString);
Tree.FullExpand;
Tree.Selected:=Tree.Items[0];
//Tree.SetFocus;
end;

{********* Pairing function :Compute the pairing according to the specified algorithm***********}
function TBNCurvePairing.Paire(P: FpPoint; Q: Fp2Point): Fp12Int;
begin
if (P.Infinity)or(Q.Infinity) then raise Exception.Create('Can''t compute pairing for infinit points ....');
case PairingAlgorithm of
  bpOptAte:Result:=OptAtePairing(P,Q);
  bpR_Ate:Result:=R_AtePairing(P,Q);
  bpEta:Result:=EtaPairing(P,Q);
  bpTate:Result:=TatePairing(P,Q);
end;

end;


{******************** Set the loop mode of the pairing :Binary / Negatve Representation*************}
procedure TBNCurvePairing.SetLoopMode(Value: LoopPoweringMode);
begin
_LoopMode:=Value;
case value of
lpmBinary:begin
          case PairingAlgorithm of
          bpOptAte:Loop:=CurveParams.LoopBin;
          bpR_Ate:Loop:=CurveParams.LoopRateBin;
          bpEta:Loop:=CurveParams.LoopEtaBin;
          bpTate:Loop:=CurveParams.LoopTateBin;
          end;

          end;
lpmNaf:begin
       case PairingAlgorithm of
       bpOptAte:Loop:=CurveParams.LoopNaf;
       bpR_Ate:Loop:=CurveParams.LoopRateNaf;
       bpEta:Loop:=CurveParams.LoopEtaNaf;
       bpTate:Loop:=CurveParams.LoopTateNaf;
       end;
       end;
lpmAuto:begin
        if CurveParams=nil then exit;
       case PairingAlgorithm of
        bpOptAte:begin
                 if ArrayHammingWeight(CurveParams.LoopNaf^)<ArrayHammingWeight(CurveParams.LoopBin^) then
                 Loop:=CurveParams.LoopNaf
                 else Loop:=CurveParams.LoopBin;
                 end;
        bpR_Ate:begin
                if ArrayHammingWeight(CurveParams.LoopRateNaf^)<ArrayHammingWeight(CurveParams.LoopRateBin^) then
                Loop:=CurveParams.LoopRateNaf
                else Loop:=CurveParams.LoopRateBin;
                end;
        bpEta:begin
              if ArrayHammingWeight(CurveParams.LoopEtaNaf^)<ArrayHammingWeight(CurveParams.LoopEtaBin^) then
              Loop:=CurveParams.LoopEtaNaf
              else Loop:=CurveParams.LoopEtaBin;
              end;
        bpTate:begin
               if ArrayHammingWeight(CurveParams.LoopTateNaf^)<ArrayHammingWeight(CurveParams.LoopTateBin^) then
               Loop:=CurveParams.LoopTateNaf
               else Loop:=CurveParams.LoopTateBin;
               end;
        end;
        end;
end;
end;



{******************** Final Exponentiation Step :Compute f^((p^12-1)/r) ***********************}
Procedure TBNCurvePairing.FinalPowerOptimalAteBN(f:Fp12Int;u:LInt;PoweringMode:FpPoweringMode;var Result:Fp12Int); // Compute f^((p^12-1)/r)
var Rst:Fp12Int;
    t:array[0..4]of Fp12Int;
begin
_Conjugate_FP12(f,t[0]);
_Inv_FP12(f,t[1]);
_Mul_FP12(t[0],t[1],Rst);
_Pow_FP12_P2(Rst,t[2]);
_Mul_FP12(t[2],Rst,Rst);
_Pow_FP12(Rst,u,t[0],PoweringMode);  ///f^x
_Sqr_FP12(t[0],t[0]);     //f^2x
_Conjugate_FP12(t[0],t[1]);  /// f-^2x
_Sqr_FP12(t[0],t[2]);       ///f^4x
_Mul_FP12(t[2],t[0],t[0]);   //// f^6x
_Pow_FP12(t[0],u ,t[2],PoweringMode); //f^6x^^2
_Sqr_FP12(t[2],t[3]);     ///f^12x^2
_Pow_FP12(t[3],u,t[4],PoweringMode);    ////f^12x3
_Mul_FP12(t[4],t[2],t[4]);
_Mul_FP12(t[4],t[0],t[4]);  //a
_Mul_FP12(t[4],t[1],t[3]);  //b
_Mul_FP12(t[4],t[2],t[1]);
_Mul_FP12(t[1],Rst,t[1]);
_Pow_FP12_P(t[3],t[0]);
_Pow_FP12_P2(t[4],t[2]);
_Conjugate_FP12(Rst,t[4]);
_Mul_FP12(t[3],t[4],t[3]);
_Pow_FP12_P3(t[3],t[3]);
_Mul_FP12(t[3],t[2],Result);
_Mul_FP12(Result,t[0],Result);
_Mul_FP12(Result,t[1],Result);
end;

{******************** Compute the Optimal Ate Pairing ***********************}
Function TBNCurvePairing.OptAtePairing(Pt:FpPoint;Qt:Fp2Point):Fp12Int;
var T,Q1,mQt,Q2,Q3,S:Fp2Point;
    i:integer;
    f:Fp12Int;
begin
f.SetTowerParams(CurveParams.TowerParam);
f.SetToOne;
T:=Qt;
_Neg_Fp2_Point(Qt,mQt);
T.SetPairingPointCoordinates(Pt.X,Pt.Y);
Q1.SetPairingPointCoordinates(Pt.X,Pt.Y);
Q2.SetPairingPointCoordinates(Pt.X,Pt.Y);
Q3.SetPairingPointCoordinates(Pt.X,Pt.Y);
S.SetPairingPointCoordinates(Pt.X,Pt.Y);
T.ComputeLigneValue:=true;
S.ComputeLigneValue:=true;
for i:=Length(Loop^)-2 Downto 0 do begin
                        case CoordinatesSystem of
                          csAffine:_Double_Affine_Fp2_Point(T,T,true);
                          csJacobian:_Double_Jacobian_Fp2_Point(T,T,True);
                          csProjective:_Double_Projective_Fp2_Point(T,T,True);
                        end;
                        _Sqr_FP12(f,f);
                        _Sparse_Mul_FP12(f,T.LineAtP,f,TwistMode);
                        if  Loop^[i]=1 then begin
                                         case CoordinatesSystem of
                                         csAffine:_Add_Affine_Fp2_Point(Qt,T,T,True);
                                         csJacobian:_Add_Jacobian_Fp2_Point(Qt,T,T,True);
                                         csProjective:_Add_Projective_Fp2_Point(Qt,T,T,True);
                                         end;
                                         _Sparse_Mul_FP12(f,T.LineAtP,f,TwistMode);
                                         end
                        else if (Loop^[i]=-1) then begin
                                                case CoordinatesSystem of
                                                csAffine:_Add_Affine_Fp2_Point(mQt,T,T,True);
                                                csJacobian:_Add_Jacobian_Fp2_Point(mQt,T,T,True);
                                                csProjective:_Add_Projective_Fp2_Point(mQt,T,T,True);
                                                end;
                                                _Sparse_Mul_FP12(f,T.LineAtP,f,TwistMode);
                                                end;
                        end;
if _IsNeg(CurveParams.u) then begin
                           _Conjugate_FP12(f,f);
                           _Neg_Fp2_Point(T,T);
                           end;
_Frobenius_Map(Qt,Q1);
_Frobenius_Map2(Qt,Q2);
_Neg_Fp2_Point(Q2,Q2);
case CoordinatesSystem of
    csAffine:_Add_Affine_Fp2_Point(Q1,T,T,True);
    csJacobian:_Add_Jacobian_Fp2_Point(Q1,T,T,True);
    csProjective:_Add_Projective_Fp2_Point(Q1,T,T,True);
end;
_Sparse_Mul_FP12(f,T.LineAtP,f,TwistMode);
case CoordinatesSystem of
    csAffine:_Add_Affine_Fp2_Point(Q2,T,T,True);
    csJacobian:_Add_Jacobian_Fp2_Point(Q2,T,T,True);
    csProjective:_Add_Projective_Fp2_Point(Q2,T,T,True);
end;
_Sparse_Mul_FP12(f,T.LineAtP,f,TwistMode);

FinalPowerOptimalAteBN(f,CurveParams.u,Fp12PoweringMode,Result);

end;


{******************** Compute the R-Ate Pairing ***********************}
Function TBNCurvePairing.R_AtePairing(Pt:FpPoint;Qt:Fp2Point):Fp12Int;
var T,mQt,tmpQt:Fp2Point;
    i:integer;
    f:Fp12Int;
begin
f.SetTowerParams(CurveParams.TowerParam);
f.SetToOne;
T:=Qt;
_Neg_Fp2_Point(Qt,mQt);
T.SetPairingPointCoordinates(Pt.X,Pt.Y);
T.ComputeLigneValue:=true;
for i:=Length(Loop^)-2 Downto 0 do begin
                        case CoordinatesSystem of
                          csAffine:_Double_Affine_Fp2_Point(T,T,true);
                          csJacobian:_Double_Jacobian_Fp2_Point(T,T,True);
                          csProjective:_Double_Projective_Fp2_Point(T,T,True);
                        end;
                        f:=f*f;
                        _Sparse_Mul_FP12(f,T.LineAtP,f,TwistMode);
                        if (Loop^[i]<>0) then begin
                                             if  Loop^[i]=1 then tmpQt:=Qt else tmpQt:=mQt;
                                             case CoordinatesSystem of
                                                csAffine:_Add_Affine_Fp2_Point(tmpQt,T,T,True);
                                                csJacobian:_Add_Jacobian_Fp2_Point(tmpQt,T,T,True);
                                                csProjective:_Add_Projective_Fp2_Point(tmpQt,T,T,True);
                                             end;
                                             _Sparse_Mul_FP12(f,T.LineAtP,f,TwistMode);
                                             end;
                        end;
FinalPowerOptimalAteBN(f,Qt.CurveParams.u,Fp12PoweringMode,Result);
end;

{******************** Compute the Eta Pairing ***********************}
Function TBNCurvePairing.EtaPairing(Pt:FpPoint;Qt:Fp2Point):Fp12Int;
/// the miller loop works on P instead of Q , evaluation of doubling and additions are in E(Fp) not in E(Fp2)
var T,mPt,tmpPt:FpPoint;
    i:integer;
    f:Fp12Int;
begin
f.SetTowerParams(CurveParams.TowerParam);
f.SetToOne;
T:=Pt;
T.ComputeLineAtQPFp12:=true;
T.ComputeLineAtQPFp6:=false;
T.LineAtQPFp12.SetTowerParams(CurveParams.TowerParam);
T.SetBNPairingPointCoordinates(qt.X,Qt.Y);
_Neg_Fp_Point(Pt,mPt);
for i:=Length(Loop^)-2 Downto 0 do begin
                        case CoordinatesSystem of
                          csAffine:_Double_Affine_Fp_Point(T,T,False);
                          csJacobian:_Double_Jacobian_Fp_Point(T,T,False);
                          csProjective:_Double_Projective_Fp_Point(T,T,False);
                        end;
                        f:=f*f;
                        f:=f*T.LineAtQPFp12;
                        if (Loop^[i]<>0) then begin
                                             if  Loop^[i]=1 then tmpPt:=Pt else tmpPt:=mPt;
                                             case CoordinatesSystem of
                                                csAffine:_Add_Affine_Fp_Point(tmpPt,T,T,False);
                                                csJacobian:_Add_Jacobian_Fp_Point(tmpPt,T,T,False);
                                                csProjective:_Add_Projective_Fp_Point(T,tmpPt,T,False);
                                             end;
                                             f:=f*T.LineAtQPFp12;
                                             end;
                        end;
if _IsNeg(curveparams.u) then _Conjugate_FP12(f,f);
FinalPowerOptimalAteBN(f,Qt.CurveParams.u,Fp12PoweringMode,Result);
end;

{******************** Compute the Tate Pairing ***********************}
Function TBNCurvePairing.TatePairing(Pt:FpPoint;Qt:Fp2Point):Fp12Int;
/// the miller loop works on P instead of Q , evaluation of doubling and additions are in E(Fp) not in E(Fp2)
var T,mPt,tmpPt:FpPoint;
    i:integer;
    f:Fp12Int;
begin
f.SetTowerParams(CurveParams.TowerParam);
f.SetToOne;
T:=Pt;
T.ComputeLineAtQPFp12:=true;
T.ComputeLineAtQPFp6:=false;
T.LineAtQPFp12.SetTowerParams(CurveParams.TowerParam);
T.SetBNPairingPointCoordinates(qt.X,Qt.Y);
_Neg_Fp_Point(Pt,mPt);
for i:=Length(Loop^)-2 Downto 0 do begin
                        case CoordinatesSystem of
                          csAffine:_Double_Affine_Fp_Point(T,T,False);
                          csJacobian:_Double_Jacobian_Fp_Point(T,T,False);
                          csProjective:_Double_Projective_Fp_Point(T,T,False);
                        end;
                        f:=f*f;
                        f:=f*T.LineAtQPFp12;
                        if (Loop^[i]<>0) then begin
                                             if  Loop^[i]=1 then tmpPt:=Pt else tmpPt:=mPt;
                                             case CoordinatesSystem of
                                                csAffine:_Add_Affine_Fp_Point(tmpPt,T,T,False);
                                                csJacobian:_Add_Jacobian_Fp_Point(tmpPt,T,T,False);
                                                csProjective:_Add_Projective_Fp_Point(T,tmpPt,T,False);
                                             end;
                                             f:=f*T.LineAtQPFp12;
                                             end;
                        end;
FinalPowerOptimalAteBN(f,Qt.CurveParams.u,Fp12PoweringMode,Result);
end;


procedure Register;
begin
  RegisterComponents('Pairings', [TBNCurvePairing]);
end;

end.
