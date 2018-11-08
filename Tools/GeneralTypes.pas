unit GeneralTypes;

interface

uses Vcl.ComCtrls, LargeIntegers,Fp6Arithmetic,Fp2Arithmetic,Fp3Arithmetic,Fp4Arithmetic,Fp8Arithmetic,Fp9Arithmetic, System.SysUtils, System.classes;


type

   CurvesFamily=(cfSuperSingular,cfBN,cfBLS12,cfBLS24,cfKSS18,cfKSS16,cfMNT,cfKSS36,cfBLS27);
   CoordinatesSystem=(csAffine,csProjective,csJacobian);
   LoopPoweringMode=(lpmBinary,lpmNaf,lpmAuto);
   TTwistModel=(twDType,twMType);
   FpPoweringMode=(pmNormal,pmBeuchat,pmKarbina);   // Diffrent Algorithms for exponentiation over the Cyclotpmic Group
   FpKSS16PoweringMode=(pmkssNormal);


   PtrCurveParams=^CurveParams;
   CurveParams=Record
                         {Definition of a Curve with respect the the Weistrass model   Y^2=X^3+A*X+B (mod P)}
               Family:CurvesFamily;    // Fmaily of the curve (Supersingular, BN, MNT, KSS,.....)
               SecurityLevel:String;    /// Security Level of the pairing on the curve
               A,B,   // parametres of the curve   y^2=x^3+Ax+B
               P,    //  the prime Modulo
               N,     // Total number of the points #E(Fp)
               R,    // Order of the sub-group on the curve (number of the pointes)
               H:LInt;    // Cofactor of the curve H=P div R
               FieldParam:PtrFieldParams;

                                          { parametres for BLS/BN curves }
               u,usqr,ucub,     // the paramater of generation for the BN curve
               Tr,   // Forbeinus Trace of the curve
               Lp,   // Size of the Miller's Loop
               Rtw:LInt;   /// Ordre of the Twist Curve
               Htw:LInt;  // Cofactor of the Twist Curve
               Atw,Btw:Fp2Int;  // Parametres of the Sextic Twist Curve
               AtwFp4,BtwFp4:Fp4Int;
               AtwFp3,BtwFp3:Fp3Int;
               AtwFp6,BtwFp6:Fp6Int;
               AtwFp9,BtwFp9:Fp9Int;
               TowerParam:PtrTowerParams12;
               TowerParam2:PtrTowerParams24;
               TowerParam3:PtrTowerParams18;
               FrobeniusMapConstX,FrobeniusMapConstY:array[0..2] of Fp2Int;// for BN/BLS12 Curves
               FrobeniusMapConstX_Fp3,FrobeniusMapConstY_Fp3:Fp3Int;// for KSS18 Curves
               FrobeniusMapConstX_Fp4,FrobeniusMapConstY_Fp4:Fp4Int;// for BLS24 Curves
               FrobeniusMapConstX_Fp6,FrobeniusMapConstY_Fp6:Fp6Int;// for KSS36 Curves
               FrobeniusMapConstX_Fp9,FrobeniusMapConstY_Fp9:Fp9Int;// for BLS27 Curves
               CoordSys:CoordinatesSystem;
               LoopBin,LoopNaf,
               LoopRateBin,LoopRateNaf,
               LoopEtaBin,LoopEtaNaf,
               LoopTateNaf,LoopTateBin:PLIntArrayForm;
               TwistMode:TTwistModel;
               BasePointX,BasePointY:Lint;
               TwistBasePointX,TwistBasePointy:Fp2Int;
               TwistBasePointX_Fp3,TwistBasePointy_Fp3:Fp3Int;
               BLS24TwistGeneratorX,BLS24TwistGeneratorY:Fp4Int;
               KSS36TwistGeneratorX,KSS36TwistGeneratorY:Fp6Int;
               BLS27TwistGeneratorX,BLS27TwistGeneratorY:Fp9Int;
               TwsistGeneratorSeed:Word;
               power:Lint;
            end;
   TCurve=class(Tcomponent)
           public
           CurveParams:PtrCurveParams;
           procedure GenerateParamsTreeView(Tree:TTreeView);virtual;
           constructor Create(AOwner : TComponent); override;
           destructor destroy;
          end;

implementation


{ TCurve }

{ TCurve }


constructor TCurve.Create(AOwner : TComponent);
begin
inherited Create(AOwner);
//initialize;
end;

destructor TCurve.destroy;
begin
inherited destroy;
end;

procedure TCurve.GenerateParamsTreeView(Tree: TTreeView);
begin

end;


end.
