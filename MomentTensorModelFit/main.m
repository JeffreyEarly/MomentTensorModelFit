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
        
        GLEquation *equation = [[GLEquation alloc] init];
        
        // This model requires initial conditions (Mxx0, Myy0, Mxy0), and a one-dimensional time variable t.
        // It takes one parameter, kappa.
        // It returns (Mxx, Myy, Mxy) for all time t.
        NSArray * (^diffusivityModel) (GLFloat, GLFloat, GLFloat, GLFunction *, GLScalar * ) = ^(GLFloat Mxx0, GLFloat Myy0, GLFloat Mxy0, GLFunction *t, GLScalar *kappa ) {
            GLFunction *Mxx = [[t times: [kappa times: @(2)]] plus: @(Mxx0)];
            GLFunction *Myy = [[t times: [kappa times: @(2)]] plus: @(Myy0)];
            GLFunction *Mxy = [[t times: @(0)] plus: @(Mxy0)];
            
            return @[Mxx, Myy, Mxy];
        };
        
		
		// This model requires initial conditions (Mxx0, Myy0, Mxy0), and a one-dimensional time variable t.
        // It takes three parameters, kappa, sigma, and theta.
        // It returns (Mxx, Myy, Mxy) for all time t.
        NSArray * (^strainDiffusivityModel) (GLFloat, GLFloat, GLFloat, GLFunction *, GLScalar *, GLScalar *, GLScalar * ) = ^(GLFloat Mxx0, GLFloat Myy0, GLFloat Mxy0, GLFunction *t, GLScalar *kappa, GLScalar *sigma, GLScalar *theta ) {
            GLScalar *cos_t = [theta cos];
            GLScalar *sin_t = [theta sin];
            
            GLScalar *cos2 = [cos_t times: cos_t];
            GLScalar *sin2 = [sin_t times: sin_t];
            GLScalar *cossin = [cos_t times: sin_t];
            
            GLScalar * tks = [[kappa times: @(2)] dividedBy: sigma];
            
            GLScalar * A = [[[cos2 times: @(Mxx0)] plus: [sin2 times: @(Myy0)]] plus: [[cossin times: @(2.*Mxy0)] plus: tks]];
            GLScalar * B = [[[sin2 times: @(Mxx0)] plus: [cos2 times: @(Myy0)]] minus: [[cossin times: @(2.*Mxy0)] plus: tks]];
            GLScalar * C = [[[cossin times: @(-Mxx0)] plus: [cossin times: @(Myy0)]] plus: [[cos2 minus: sin2] times: @(Mxy0)]];
            
            GLFunction *Maa = [[[[t times: sigma] exponentiate] times: A] minus: tks];
            GLFunction *Mbb = [[[[t times: [sigma negate]] exponentiate] times: B] plus: tks];
            GLFunction *Mab = [[t scalarMultiply: 0.0] plus: C];
            
            GLFunction *Mxx = [[[Maa times: cos2] plus: [Mbb times: sin2]] plus: [Mab times: [cossin times: @(-2.)]]];
            GLFunction *Myy = [[[Maa times: sin2] plus: [Mbb times: cos2]] plus: [Mab times: [cossin times: @(2.)]]];
            GLFunction *Mxy = [[[Maa times: cossin] minus: [Mbb times: cossin]] plus: [Mab times: [cos2 minus: sin2]]];
            
            return @[Mxx, Myy, Mxy];
        };
        
        // This model requires initial conditions (Mxx0, Myy0, Mxy0), and a one-dimensional time variable t.
        // It takes three parameters, kappa, sigma, and theta.
        // It returns (Mxx, Myy, Mxy) for all time t.
        NSArray * (^strainVorticityDiffusivityModel) (GLFloat, GLFloat, GLFloat, GLFunction *, GLScalar *, GLScalar *, GLScalar *, GLScalar * ) = ^(GLFloat Mxx0, GLFloat Myy0, GLFloat Mxy0, GLFunction *t, GLScalar *kappa, GLScalar *sigma, GLScalar *theta, GLScalar *zeta ) {
            GLScalar *cos_t = [theta cos];
            GLScalar *sin_t = [theta sin];
            
            GLScalar *cos2 = [cos_t times: cos_t];
            GLScalar *sin2 = [sin_t times: sin_t];
            GLScalar *cossin = [cos_t times: sin_t];
            
            GLScalar * s2 = [[sigma times: sigma] minus: [zeta times: zeta]];
            GLScalar * s = [s2 sqrt];
            GLScalar * tks = [[[kappa times: @(2)] times: sigma] dividedBy: s2];
            GLScalar * sig_s = [sigma dividedBy: s];
            GLScalar * zeta_s = [zeta dividedBy: s];
            
            GLScalar * one_plus_sigma_over_s_div2 = [[sig_s plus: @(1)] times: @(0.5)];
            GLScalar * one_minus_sigma_over_s_div2 = [[sig_s minus: @(1)] times: @(-0.5)];
            GLFunction * exp_s_t = [[t times: s] exponentiate];
            GLFunction * exp_minus_s_t = [[t times: [s negate]] exponentiate];
            
            GLScalar * Mxx1 = [[[cos2 times: @(Mxx0)] plus: [sin2 times: @(Myy0)]] plus: [cossin times: @(2.*Mxy0)]];
            GLScalar * Myy1 = [[[sin2 times: @(Mxx0)] plus: [cos2 times: @(Myy0)]] minus: [cossin times: @(2.*Mxy0)]];
            GLScalar * Mxy1 = [[[cossin times: @(-Mxx0)] plus: [cossin times: @(Myy0)]] plus: [[cos2 minus: sin2] times: @(Mxy0)]];
            
            GLScalar * A = [[[[one_plus_sigma_over_s_div2 times: Mxx1] minus: [zeta_s times: Mxy1]] minus: [one_minus_sigma_over_s_div2 times: Myy1]] plus: tks];
            GLScalar * B = [[[[one_minus_sigma_over_s_div2 times: Mxx1] plus: [zeta_s times: Mxy1]] minus: [one_plus_sigma_over_s_div2 times: Myy1]] plus: tks];
            GLScalar * C = [[[zeta_s negate] times: [Mxx1 plus: Myy1]] plus: [[sig_s times: @2] times: Mxy1]];
            
            GLFunction *Maa = [[[exp_s_t times: one_plus_sigma_over_s_div2] times: A] plus: [[exp_minus_s_t times: one_minus_sigma_over_s_div2] times: B]];
            Maa = [Maa plus: [[C times: zeta_s] times: @(0.5)]];
            Maa = [Maa minus: [[[[zeta times: zeta] times: t] plus: sigma] times: [[kappa times: @2] dividedBy: s2]]];
            GLFunction *Mbb = [[[[exp_s_t times: one_minus_sigma_over_s_div2] times: A] plus: [[exp_minus_s_t times: one_plus_sigma_over_s_div2] times: B]] negate];
            Mbb = [Mbb plus: [[C times: zeta_s] times: @(0.5)]];
            Mbb = [Mbb minus: [[[[zeta times: zeta] times: t] minus: sigma] times: [[kappa times: @2] dividedBy: s2]]];
            GLFunction *Mab = [[[[exp_s_t times: zeta_s] times: A] minus: [[exp_minus_s_t times: zeta_s] times: B]] times: @(0.5)];
            Mab = [Mab plus: [[C times: @(0.5)] times: sig_s]];
            Mab = [Mab minus: [[[zeta times: sigma] times: t] times: [[kappa times: @2] dividedBy: s2]]];
            
            GLFunction *Mxx = [[[Maa times: cos2] plus: [Mbb times: sin2]] plus: [Mab times: [cossin times: @(-2.)]]];
            GLFunction *Myy = [[[Maa times: sin2] plus: [Mbb times: cos2]] plus: [Mab times: [cossin times: @(2.)]]];
            GLFunction *Mxy = [[[Maa times: cossin] minus: [Mbb times: cossin]] plus: [Mab times: [cos2 minus: sin2]]];
            
            return @[Mxx, Myy, Mxy];
        };
        
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
        
        if (0) {
            GLDimension *tDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 10 domainMin: 0 length: 9*1800];
            GLFunction *t = [GLFunction functionOfRealTypeFromDimension: tDim withDimensions: @[tDim] forEquation: equation];
            GLScalar *kappa = [GLScalar scalarWithValue: 0.2014 forEquation: equation];
            GLScalar *sigma = [GLScalar scalarWithValue: 3.4342e-06 forEquation: equation];
            GLScalar *theta = [GLScalar scalarWithValue: -32.2493*M_PI/180 forEquation: equation];
            GLScalar *zeta = [GLScalar scalarWithValue: -3.0000e-06 forEquation: equation];
            //GLScalar *zeta = [GLScalar scalarWithValue: 1e-9 forEquation: equation];
            GLFloat Mxx0 = 1.8499e+05;
            GLFloat Myy0 = 2.7761e+05;
            GLFloat Mxy0 = -5.0412e+04;
            
            NSArray *tensorComps = strainVorticityDiffusivityModel( Mxx0, Myy0, Mxy0, t, kappa, sigma, theta, zeta);
            //NSArray *tensorComps = strainDiffusivityModel( Mxx0, Myy0, Mxy0, t, kappa, sigma, theta);
            GLFunction *Mxx=tensorComps[0]; GLFunction *Myy=tensorComps[1]; GLFunction *Mxy=tensorComps[2];
            for (int i=0; i<tDim.nPoints; i++) {
                printf("(Mxx,Myy,Mxy)=(%f, %f, %f)\n", Mxx.pointerValue[i], Myy.pointerValue[i], Mxy.pointerValue[i]);
            }
            printf("\n\n");
            NSArray *ellipseComps = tensorCompsToEllipseComps( tensorComps );
            GLFunction *semiMajor=ellipseComps[0]; GLFunction *semiMinor=ellipseComps[1]; GLFunction *angle=ellipseComps[2];
            for (int i=0; i<tDim.nPoints; i++) {
                printf("(a,b,theta)=(%f, %f, %f)\n", semiMajor.pointerValue[i], semiMinor.pointerValue[i], angle.pointerValue[i]*180/M_PI);
            }
            
            return 0;
        }
        
		
		NSFileManager *fileManager = [[NSFileManager alloc] init];
		NSString *folderPath = @"/Users/jearly/Documents/LatMix/drifters/ObservationalData/griddedRhoDrifterMomementEllipses/";
        folderPath = @"/Users/jearly/Documents/LatMix/drifters/synthetic/moment-ellipses/synthetic-diffusive/";
        folderPath = @"/Users/jearly/Documents/LatMix/drifters/synthetic/moment-ellipses/synthetic-strained-diffusive/";
		folderPath = @"/Users/jearly/Documents/LatMix/drifters/synthetic/moment-ellipses/synthetic-vorticity-strained-diffusive/";
		NSArray *ellipseFiles = [fileManager contentsOfDirectoryAtPath: folderPath error: nil];
		
		NSMutableString *outputData = [NSMutableString stringWithFormat: @""];
		NSUInteger i=1;
		
		for ( NSString *filename in ellipseFiles)
		{
			if ( ![filename.pathExtension isEqualToString: @"txt"]) continue;
			MomentEllipseFileReader *file = [[MomentEllipseFileReader alloc] initWithURL: [NSURL fileURLWithPath: [NSString stringWithFormat: @"%@%@", folderPath, filename]] equation: equation];
			if (!file) continue;
			
			NSArray *drifterIDs = [[filename stringByDeletingPathExtension] componentsSeparatedByString: @"_"];
			NSMutableString *idString = [NSMutableString stringWithFormat: @"drifterIDs{%lu}=[", i];
			for (NSString *name in drifterIDs) {
				if ([drifterIDs indexOfObject: name] == 0) {
					continue; // the first string is just 'drifter'.
				}
				else if ([drifterIDs indexOfObject: name] == drifterIDs.count-1) {
					[idString appendFormat: @"%@", name];
				} else {
					[idString appendFormat: @"%@;", name];
				}
			}
			[idString appendString: @"]; "];
			[outputData appendString: idString];
			
			
			GLFloat kappaScale = 0.1;
			GLScalar *kappa = [GLScalar scalarWithValue: log(0.1/kappaScale) forEquation: equation];
			GLScalar *kappaDelta = [GLScalar scalarWithValue: 2.0 forEquation: equation];
			
			GLMinimizationOperation *minimizer = [[GLMinimizationOperation alloc] initAtPoint: @[kappa] withDeltas: @[kappaDelta] forFunction:^(NSArray *xArray) {
				GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
				NSArray *tensorComps = diffusivityModel( file.Mxx0, file.Myy0, file.Mxy0, file.t, kappaUnscaled);
				NSArray *ellipseComps = tensorCompsToEllipseComps( tensorComps );
		 
				EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[file.a, file.b, file.angle]];
			
				return error.result[0];
			}];
			
			GLScalar *minKappa = [[minimizer.result[0] exponentiate] times: @(kappaScale)];
			GLScalar *minError = minimizer.result[1];
			
			[outputData appendFormat: @"model1_error(%lu)=%g; model1_kappa(%lu)=%g; ", i, *(minError.pointerValue), i, *(minKappa.pointerValue)];
			
			NSLog(@"%@---diffusivity model total error: %f @ (kappa)=(%.4f)", filename.lastPathComponent, *(minError.pointerValue), *(minKappa.pointerValue));
			
            
            
            
			// kappaDelta is carefully chosen. It needs to represent a sort of 'size of parameter space' that we want to explore.
			// So here, we make it move around in fairly big chunks, like orders of magnitude.
			// Yup, big steps are the best, for ALL parameters.
			kappa = [GLScalar scalarWithValue: log(.1/kappaScale) forEquation: equation];
			kappaDelta = [GLScalar scalarWithValue: 3.0 forEquation: equation];
			
			GLFloat sigmaScale = 4.0E-6;
			GLScalar *sigma = [GLScalar scalarWithValue: log(1E-6/sigmaScale) forEquation: equation];
			GLScalar *sigmaDelta = [GLScalar scalarWithValue: 0.5 forEquation: equation];
			
			GLScalar *theta = [GLScalar scalarWithValue: 0.0*M_PI/180. forEquation: equation];
			GLScalar *thetaDelta = [GLScalar scalarWithValue: 45.*M_PI/180. forEquation: equation];
			
			minimizer = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, sigma, theta] withDeltas: @[kappaDelta, sigmaDelta, thetaDelta] forFunction:^(NSArray *xArray) {
				GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
				GLScalar *sigmaUnscaled = [[xArray[1] exponentiate] times: @(sigmaScale)];
				NSArray *tensorComps = strainDiffusivityModel( file.Mxx0, file.Myy0, file.Mxy0, file.t, kappaUnscaled, sigmaUnscaled, xArray[2]);
				NSArray *ellipseComps = tensorCompsToEllipseComps( tensorComps );
				
				EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[file.a, file.b, file.angle]];
				
				return error.result[0];
			}];
			
            NSArray *results = minimizer.result;
            
			minKappa = [[minimizer.result[0] exponentiate] times: @(kappaScale)];
			GLScalar *minSigma = [[minimizer.result[1] exponentiate] times: @(sigmaScale)];
			GLScalar *minTheta = minimizer.result[2];
			minError = minimizer.result[3];
			
			[outputData appendFormat: @"model2_error(%lu)=%g; model2_kappa(%lu)=%g; model2_sigma(%lu)=%g; model2_theta(%lu)=%g;\n", i, *(minError.pointerValue), i, *(minKappa.pointerValue), i, *(minSigma.pointerValue), i, *(minTheta.pointerValue)];
			
			NSLog(@"%@---strain-diffusivity model total error: %f (kappa,sigma,theta)=(%.4f,%.3g,%.1f)", filename.lastPathComponent, *(minError.pointerValue), *(minKappa.pointerValue),*(minSigma.pointerValue),(*(minTheta.pointerValue))*180./M_PI);
            
            
            sigmaScale = 1e-5;
            sigma = [GLScalar scalarWithValue: log(1E-5/sigmaScale) forEquation: equation];
            
			
			
            GLFloat sScale = 1E-6;
			GLScalar *s = [GLScalar scalarWithValue: log(1e-6/sScale) forEquation: equation];
			GLScalar *sDelta = [GLScalar scalarWithValue: 0.5 forEquation: equation];
			
			GLScalar *alpha = [GLScalar scalarWithValue: 0 forEquation: equation];
			GLScalar *alphaDelta = [GLScalar scalarWithValue: 1.0 forEquation: equation];
            
            // initializing with the results from the previous calculation. we can use log searches if we look for both positive and negative vorticity.
            GLMinimizationOperation *minimizer_positive = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, s, theta, alpha] withDeltas: @[kappaDelta, sDelta, thetaDelta, alphaDelta] forFunction:^(NSArray *xArray) {
				GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
				GLScalar *sUnscaled = [[xArray[1] exponentiate] times: @(sScale)];
				GLScalar *sigmaUnscaled = [sUnscaled times: [xArray[3] cosh]];
                GLScalar *zetaUnscaled = [sUnscaled times: [xArray[3] sinh]];
				NSArray *tensorComps = strainVorticityDiffusivityModel( file.Mxx0, file.Myy0, file.Mxy0, file.t, kappaUnscaled, sigmaUnscaled, xArray[2], zetaUnscaled);
				NSArray *ellipseComps = tensorCompsToEllipseComps( tensorComps );
				
				EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[file.a, file.b, file.angle]];
				
				return error.result[0];
			}];
            NSArray *posResults = minimizer_positive.result;
            
			alphaDelta = [GLScalar scalarWithValue: -1.0 forEquation: equation];
			GLMinimizationOperation *minimizer_negative = [[GLMinimizationOperation alloc] initAtPoint: @[kappa, s, theta, alpha] withDeltas: @[kappaDelta, sDelta, thetaDelta, alphaDelta] forFunction:^(NSArray *xArray) {
				GLScalar *kappaUnscaled = [[xArray[0] exponentiate] times: @(kappaScale)];
				GLScalar *sUnscaled = [[xArray[1] exponentiate] times: @(sScale)];
				GLScalar *sigmaUnscaled = [sUnscaled times: [xArray[3] cosh]];
                GLScalar *zetaUnscaled = [sUnscaled times: [xArray[3] sinh]];
				NSArray *tensorComps = strainVorticityDiffusivityModel( file.Mxx0, file.Myy0, file.Mxy0, file.t, kappaUnscaled, sigmaUnscaled, xArray[2], zetaUnscaled);
				NSArray *ellipseComps = tensorCompsToEllipseComps( tensorComps );
				
				EllipseErrorOperation *error = [[EllipseErrorOperation alloc] initWithParametersFromEllipseA: ellipseComps ellipseB:@[file.a, file.b, file.angle]];
				
				return error.result[0];
			}];
            NSArray *negResults = minimizer_negative.result;
            
			GLScalar *minPos = posResults[4];
			GLScalar *minNeg = negResults[4];
			GLFloat pos = *(minPos.pointerValue);
			GLFloat neg = *(minNeg.pointerValue);
			
			results = pos < neg ? posResults : negResults;
			
			minKappa = [[results[0] exponentiate] times: @(kappaScale)];
			GLScalar *minS = [[results[1] exponentiate] times: @(sScale)];
			minSigma = [minS times: [results[3] cosh]];
            minTheta = results[2];
            GLScalar *minZeta = [minS times: [results[3] sinh]];
			minError = results[4];
			
			[outputData appendFormat: @"model3_error(%lu)=%g; model3_kappa(%lu)=%g; model3_sigma(%lu)=%g; model3_theta(%lu)=%g; model3_zeta(%lu)=%g;\n", i, *(minError.pointerValue), i, *(minKappa.pointerValue), i, *(minSigma.pointerValue), i, *(minTheta.pointerValue), i, *(minZeta.pointerValue)];
			
			NSLog(@"%@---strain-vorticity-diffusivity model total error: %f (kappa,sigma,theta)=(%.4f,%.3g,%.1f,%.3g)", filename.lastPathComponent, *(minError.pointerValue), *(minKappa.pointerValue),*(minSigma.pointerValue),(*(minTheta.pointerValue))*180./M_PI, *(minZeta.pointerValue));
            
			i++;
        }
		
		[outputData writeToFile: [NSString stringWithFormat: @"%@BestFitParameters_area_divergence.m", folderPath] atomically: YES encoding: NSUTF8StringEncoding error: nil];
    }
    return 0;
}
