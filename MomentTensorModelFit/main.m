//
//  main.m
//  MomentTensorFit
//
//  Created by Jeffrey Early on 9/16/13.
//
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import <GLNumericalModelingKit/GLOperationOptimizer.h>
#import "MomentEllipseFileReader.h"
#import "EllipseErrorOperation.h"

int main(int argc, const char * argv[])
{
	
    @autoreleasepool {
        
        NSURL *ellipseParameterFile = [NSURL fileURLWithPath:@"/Users/jearly/Dropbox/Documents/Projects/LatMix/umassd_drifters/ObservationalData/griddedRho1Drifters_moment_ellipse.txt"];
        GLEquation *equation = [[GLEquation alloc] init];
        
        MomentEllipseFileReader *file = [[MomentEllipseFileReader alloc] initWithURL: ellipseParameterFile equation: equation];
        
//        GLFloat sigma = 3.77732e-06;
//        GLFloat theta = -32.4243*M_PI/180.;
//        GLFloat kappa = 0.238886;
        
        
        // This model requires initial conditions (Mxx0, Myy0, Mxy0), and a one-dimensional time variable t.
        // It takes one parameter, kappa.
        // It returns (Mxx, Myy, Mxy) for all time t.
        NSArray * (^diffusivityModel) (GLFloat, GLFloat, GLFloat, GLFunction *, GLScalar * ) = ^(GLFloat Mxx0, GLFloat Myy0, GLFloat Mxy0, GLFunction *t, GLScalar *kappa ) {
            GLFunction *Mxx = [[t times: [kappa times: @(2)]] plus: @(Mxx0)];
            GLFunction *Myy = [[t times: [kappa times: @(2)]] plus: @(Myy0)];
            GLFunction *Mxy = [[t times: @(0)] plus: @(Mxy0)];
            
            return @[Mxx, Myy, Mxy];
        };
        
//        NSArray * (^strainDiffusivityModel) (GLFloat, GLFloat, GLFloat, GLFunction *, GLScalar *, GLScalar *, GLScalar * ) = ^(GLFloat Mxx0, GLFloat Myy0, GLFloat Mxy0, GLFunction *t, GLScalar *kappa, GLScalar *sigma, GLScalar *theta ) {
//            GLFunction *cos_t = [theta cos];
//            GLFunction *sin_t = [theta sin];
//            
//            GLFunction * cos2 = [cos_t times: cos_t];
//            GLFunction * sin2 = [sin_t times: sin_t];
//            GLFunction * cossin = [cos_t times: sin_t];
//            
//            GLFunction * tks = [[kappa times: @(2)] dividedBy: sigma];
//            
//            GLFunction * A = [[[cos2 scalarMultiply: Mxx0] plus: [sin2 scalarMultiply: Myy0]] plus: [[cossin scalarMultiply: 2.*Mxy0] plus: tks]];
//            GLFunction * B = [[[sin2 scalarMultiply: Mxx0] plus: [cos2 scalarMultiply: Myy0]] minus: [[cossin scalarMultiply: 2.*Mxy0] plus: tks]];
//            GLFunction * C = [[[cossin scalarMultiply: -Mxx0] plus: [cossin scalarMultiply: Myy0]] plus: [[cos2 minus: sin2] scalarMultiply: Mxy0]];
//            
//            GLFunction *Maa = [[[[t times: sigma] exponentiate] times: A] minus: tks];
//            GLFunction *Mbb = [[[[t times: [sigma negate]] exponentiate] times: B] plus: tks];
//            GLFunction *Mab = [[t scalarMultiply: 0.0] plus: C];
//            
//            GLFunction *Mxx = [[[Maa times: cos2] plus: [Mbb times: sin2]] plus: [Mab times: [cossin scalarMultiply: -2.]]];
//            GLFunction *Myy = [[[Maa times: sin2] plus: [Mbb times: cos2]] plus: [Mab times: [cossin scalarMultiply: 2.]]];
//            GLFunction *Mxy = [[[Maa times: cossin] minus: [Mbb times: cossin]] plus: [Mab times: [cos2 minus: sin2]]];
//            
//            return @[Mxx, Myy, Mxy];
//        };
        
        // Simple little utility function that converts tensor components, into ellipse components
        NSArray * (^tensorCompsToEllipseComps) (NSArray *) = ^( NSArray *tensorComp ) {
            GLFunction *Mxx = tensorComp[0];
            GLFunction *Myy = tensorComp[1];
            GLFunction *Mxy = tensorComp[2];
            
            GLFunction *descriminant = [[[[Mxx minus: Myy] times: [Mxx minus: Myy]] plus: [[Mxy scalarMultiply: 2] times:[Mxy scalarMultiply: 2]]] sqrt];
            GLFunction *D2_plus = [[[Mxx plus: Myy] plus: descriminant] scalarMultiply: 0.5];
            GLFunction *D2_minus = [[[Mxx plus: Myy] minus: descriminant] scalarMultiply: 0.5];
            
            GLFunction *semiMajor = [D2_plus sqrt];
            GLFunction *semiMinor = [D2_minus sqrt];
            GLFunction *angle = [[[D2_plus minus: Mxx] dividedBy: Mxy] atan];
            
            return @[semiMajor, semiMinor, angle];
        };
        
        GLScalar *kappa = [GLScalar scalarWithValue: 0.1 forEquation: equation];
        GLScalar *kappaDelta = [GLScalar scalarWithValue: 0.025 forEquation: equation];
        
        GLMinimizationOperation *minimizer = [[GLMinimizationOperation alloc] initAtPoint: @[kappa] withDeltas: @[kappaDelta] forFunction:^(NSArray *xArray) {
            NSArray *tensorComps = diffusivityModel( file.Mxx0, file.Myy0, file.Mxy0, file.t, xArray[0]);
            NSArray *ellipseComps = tensorCompsToEllipseComps( tensorComps );
     
            EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[file.a, file.b, file.angle]];
        
            return error.result[0];
        }];
        
        GLScalar *minPoint = minimizer.result[0];
        GLScalar *minValue = minimizer.result[1];
        
        NSLog(@"total error: %f @ (kappa)=(%.4f)", *(minValue.pointerValue), *(minPoint.pointerValue));
        
        //NSLog(@"total error: %f (kappa,sigma,theta)=(%.4f,%.3g,%.1f)", minError, bestKappa,bestSigma,bestTheta*180./M_PI);
        
    }
    return 0;
}
