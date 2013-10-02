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
#include "ellipse_ellipse_overlap.h"

int main(int argc, const char * argv[])
{
	
    @autoreleasepool {
        
        NSURL *ellipseParameterFile = [NSURL fileURLWithPath:@"/Users/jearly/Dropbox/Documents/Projects/LatMix/umassd_drifters/ObservationalData/griddedRho1Drifters_moment_ellipse.txt"];
        NSError *anError;
        NSString *fileContents = [NSString stringWithContentsOfURL: ellipseParameterFile encoding: NSASCIIStringEncoding error: &anError];
        NSUInteger numberOfLines, index, stringLength = [fileContents length];
        for (index = 0, numberOfLines = 0; index < stringLength; numberOfLines++)
            index = NSMaxRange([fileContents lineRangeForRange:NSMakeRange(index, 0)]);
		
        GLEquation *equation = [[GLEquation alloc] init];
        GLDimension *tDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: numberOfLines domainMin: 0 length: numberOfLines];
        
        GLVariable *t = [GLVariable variableOfRealTypeWithDimensions: @[tDim] forEquation:equation];
        GLVariable *a_obs = [GLVariable variableOfRealTypeWithDimensions: @[tDim] forEquation:equation];
        GLVariable *b_obs = [GLVariable variableOfRealTypeWithDimensions: @[tDim] forEquation:equation];
        GLVariable *angle_obs = [GLVariable variableOfRealTypeWithDimensions: @[tDim] forEquation:equation];
        
        NSUInteger i=0;
        NSScanner *theScanner = [NSScanner scannerWithString:fileContents];
        while (theScanner.isAtEnd == NO) {
            NSDecimal val;
            [theScanner scanDecimal:&val];
            NSDecimalNumber *number = [NSDecimalNumber decimalNumberWithDecimal: val];
            t.pointerValue[i] = number.doubleValue;
            [theScanner scanDecimal:&val];
            number = [NSDecimalNumber decimalNumberWithDecimal: val];
            a_obs.pointerValue[i] = number.doubleValue;
            [theScanner scanDecimal:&val];
            number = [NSDecimalNumber decimalNumberWithDecimal: val];
            b_obs.pointerValue[i] = number.doubleValue;
            [theScanner scanDecimal:&val];
            number = [NSDecimalNumber decimalNumberWithDecimal: val];
            angle_obs.pointerValue[i] = number.doubleValue;
            i++;
        }
        
        GLFloat D2p = pow(a_obs.pointerValue[0],2.0);
        GLFloat D2m = pow(b_obs.pointerValue[0],2.0);
        GLFloat alpha = angle_obs.pointerValue[0];
        GLFloat Mxx0 = D2p*cos(alpha)*cos(alpha) + D2m*sin(alpha)*sin(alpha);
        GLFloat Myy0 = D2p*sin(alpha)*sin(alpha) + D2m*cos(alpha)*cos(alpha);
        GLFloat Mxy0 = (D2p-D2m)*sin(alpha)*cos(alpha);
        
		//        GLFloat sigma = 3.77732e-06;
		//        GLFloat theta = -32.4243*M_PI/180.;
		//        GLFloat kappa = 0.238886;
        
        GLVariable *sigma = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: equation];
        GLVariable *theta = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: equation];
        GLVariable *kappa = [GLVariable variableOfRealTypeWithDimensions: @[] forEquation: equation];
        
        sigma.pointerValue[0] = 3.77732e-06;
        theta.pointerValue[0] = -32.4243*M_PI/180.;
        kappa.pointerValue[0] = 0.238886;
        
        GLVariable *Mxx;
        GLVariable *Myy;
        GLVariable *Mxy;
        if (0.0) { // sigma == 0.0
            Mxx = [[t times: [kappa scalarMultiply: 2]] scalarAdd: Mxx0];
            Myy = [[t times: [kappa scalarMultiply: 2]] scalarAdd: Myy0];
            Mxy = [[t scalarMultiply: 0] scalarAdd: Mxy0];
        } else {
            GLVariable *cos_t = [theta cos];
            GLVariable *sin_t = [theta sin];
            
            GLVariable * cos2 = [cos_t times: cos_t];
            GLVariable * sin2 = [sin_t times: sin_t];
            GLVariable * cossin = [cos_t times: sin_t];
            
            GLVariable * tks = [[kappa scalarMultiply: 2] dividedBy: sigma];
            
            GLVariable * A = [[[cos2 scalarMultiply: Mxx0] plus: [sin2 scalarMultiply: Myy0]] plus: [[cossin scalarMultiply: 2.*Mxy0] plus: tks]];
            GLVariable * B = [[[sin2 scalarMultiply: Mxx0] plus: [cos2 scalarMultiply: Myy0]] minus: [[cossin scalarMultiply: 2.*Mxy0] plus: tks]];
            GLVariable * C = [[[cossin scalarMultiply: -Mxx0] plus: [cossin scalarMultiply: Myy0]] plus: [[cos2 minus: sin2] scalarMultiply: Mxy0]];
            
            GLVariable *Maa = [[[[t times: sigma] exponentiate] times: A] minus: tks];
            GLVariable *Mbb = [[[[t times: [sigma negate]] exponentiate] times: B] plus: tks];
            GLVariable *Mab = [[t scalarMultiply: 0.0] plus: C];
            
            Mxx = [[[Maa times: cos2] plus: [Mbb times: sin2]] plus: [Mab times: [cossin scalarMultiply: -2.]]];
            Myy = [[[Maa times: sin2] plus: [Mbb times: cos2]] plus: [Mab times: [cossin scalarMultiply: 2.]]];
            Mxy = [[[Maa times: cossin] minus: [Mbb times: cossin]] plus: [Mab times: [cos2 minus: sin2]]];
            
        }
        
        GLVariable *descriminant = [[[[Mxx minus: Myy] times: [Mxx minus: Myy]] plus: [[Mxy scalarMultiply: 2] times:[Mxy scalarMultiply: 2]]] sqrt];
        GLVariable *D2_plus = [[[Mxx plus: Myy] plus: descriminant] scalarMultiply: 0.5];
        GLVariable *D2_minus = [[[Mxx plus: Myy] minus: descriminant] scalarMultiply: 0.5];
        
        GLVariable *semiMajor = [D2_plus sqrt];
        GLVariable *semiMinor = [D2_minus sqrt];
        GLVariable *angle = [[[D2_plus minus: Mxx] dividedBy: Mxy] atan];
        
        NSArray *topVariables = @[sigma, theta, kappa];
        NSArray *bottomVariables = @[semiMajor, semiMinor, angle];
        GLOperationOptimizer *optimizer = [[GLOperationOptimizer alloc] initWithTopVariables: topVariables bottomVariables: bottomVariables];
        unaryVectorOperation operation = optimizer.unaryVectorOperationBlock;
        NSArray *inputBuffer = @[sigma.data, theta.data, kappa.data];
        NSArray *resultBuffer = @[semiMajor.data, semiMinor.data, angle.data];
        
        kappa.pointerValue[0] = 0.1;
        operation(resultBuffer, inputBuffer);
        
        double totalObsArea = 0.0;
        for (NSUInteger iTime = 0; iTime<tDim.nPoints; iTime++) {
            double obsArea =  M_PI * a_obs.pointerValue[iTime] *  b_obs.pointerValue[iTime];
            totalObsArea += obsArea;
        }
        
        double minError = 1e10;
        double bestKappa, bestSigma, bestTheta;
        GLFloat *kappaPtr = kappa.pointerValue;
        GLFloat *sigmaPtr = sigma.pointerValue;
        GLFloat *thetaPtr = theta.pointerValue;
        GLFloat *a_obs_ptr = a_obs.pointerValue;
        GLFloat *b_obs_ptr = b_obs.pointerValue;
        GLFloat *angle_obs_ptr = angle_obs.pointerValue;
        GLFloat *a_model_ptr = semiMajor.pointerValue;
        GLFloat *b_model_ptr = semiMinor.pointerValue;
        GLFloat *angle_model_ptr = angle.pointerValue;
        
		GLVariable *errorAtTime = [GLVariable variableOfRealTypeWithDimensions: @[tDim] forEquation:equation];
		[errorAtTime zero];
		GLFloat *errorAtTimePtr = errorAtTime.pointerValue;
		
        for (kappaPtr[0] = 0.045; kappaPtr[0] < 0.055; kappaPtr[0] += 0.00025) {
            for (sigmaPtr[0] = 4e-6; sigmaPtr[0] < 5e-6; sigmaPtr[0] += 0.025e-6) {
                for (thetaPtr[0] = -35.*M_PI/180.; thetaPtr[0] < -30.*M_PI/180.; thetaPtr[0] += .1*M_PI/180.) {
					
					
					
					operation(resultBuffer, inputBuffer);
                    
					dispatch_apply(tDim.nPoints, dispatch_get_global_queue (DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^(size_t iTime) {
                        double X[4], Y[4];
                        int nroots, rtnCode;
                        double area = ellipse_ellipse_overlap(angle_obs_ptr[iTime], a_obs_ptr[iTime], b_obs_ptr[iTime], 0.0, 0.0, angle_model_ptr[iTime], a_model_ptr[iTime], b_model_ptr[iTime],0.0, 0.0,X, Y, &nroots, &rtnCode);
                        double obsArea =  M_PI * a_obs_ptr[iTime] *  b_obs_ptr[iTime];
                        double modelArea = M_PI * a_model_ptr[iTime] *  b_model_ptr[iTime];
                        errorAtTimePtr[iTime] =(modelArea - area) + (obsArea-area);
					});
					
                    
					double totalError = 0.0;
					for (NSUInteger iTime = 0; iTime<tDim.nPoints; iTime++) {
						totalError +=errorAtTimePtr[iTime];
					}
					
					totalError/=totalObsArea;
                    if (totalError < minError) {
                        bestKappa = kappaPtr[0];
                        bestSigma = sigmaPtr[0];
                        bestTheta = thetaPtr[0];
                        minError = totalError;
                    }
                }
            }
        }
		
        NSLog(@"total error: %f (kappa,sigma,theta)=(%.4f,%.3g,%.1f)", minError, bestKappa,bestSigma,bestTheta*180./M_PI);
        
    }
    return 0;
}
