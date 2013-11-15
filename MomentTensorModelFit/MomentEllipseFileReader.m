//
//  MomentEllipseFileReader.m
//  MomentTensorModelFit
//
//  Created by Jeffrey Early on 11/15/13.
//  Copyright (c) 2013 Jeffrey J. Early. All rights reserved.
//

#import "MomentEllipseFileReader.h"

@implementation MomentEllipseFileReader

- (id) initWithURL: (NSURL *) anURL equation: (GLEquation *) equation
{
    if ((self=[super init])) {
        self.file = anURL;
        NSError *anError;
        NSString *fileContents = [NSString stringWithContentsOfURL: self.file encoding: NSASCIIStringEncoding error: &anError];
        NSUInteger numberOfLines, index, stringLength = [fileContents length];
        for (index = 0, numberOfLines = 0; index < stringLength; numberOfLines++)
            index = NSMaxRange([fileContents lineRangeForRange:NSMakeRange(index, 0)]);
		
        GLDimension *tDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: numberOfLines domainMin: 0 length: numberOfLines];
        
        GLFunction *t = [GLFunction functionOfRealTypeWithDimensions: @[tDim] forEquation:equation];
        GLFunction *a_obs = [GLFunction functionOfRealTypeWithDimensions: @[tDim] forEquation:equation];
        GLFunction *b_obs = [GLFunction functionOfRealTypeWithDimensions: @[tDim] forEquation:equation];
        GLFunction *angle_obs = [GLFunction functionOfRealTypeWithDimensions: @[tDim] forEquation:equation];
        
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
        
        self.Mxx0 = D2p*cos(alpha)*cos(alpha) + D2m*sin(alpha)*sin(alpha);
        self.Myy0 = D2p*sin(alpha)*sin(alpha) + D2m*cos(alpha)*cos(alpha);
        self.Mxy0 = (D2p-D2m)*sin(alpha)*cos(alpha);
        
        self.t = t;
        self.a = a_obs;
        self.b = b_obs;
        self.angle = angle_obs;
    }
    return self;
}

@end
