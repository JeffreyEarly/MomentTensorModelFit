//
//  DrifterTracksFileReader.m
//  MomentTensorModelFit
//
//  Created by Jeffrey J. Early on 5/23/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import "DrifterTracksFileReader.h"

@implementation DrifterTracksFileReader

- (id) initWithURL: (NSURL *) anURL equation: (GLEquation *) equation
{
    if ((self=[super init])) {
        self.file = anURL;
		
        NSError *anError;
        NSString *fileContents = [NSString stringWithContentsOfURL: self.file encoding: NSASCIIStringEncoding error: &anError];
		if (anError) {
			NSLog(@"%@", [anError localizedDescription]);
			return nil;
		}
		
		// Count the number of lines, which gives use the number of time points.
        NSUInteger numberOfLines, index, stringLength = [fileContents length];
        for (index = 0, numberOfLines = 0; index < stringLength; numberOfLines++) {
            index = NSMaxRange([fileContents lineRangeForRange:NSMakeRange(index, 0)]);
		}
		
		// Count the number of tabs, which gives us the number of drifters
		NSScanner *theScanner = [NSScanner scannerWithString:fileContents];
		NSString *firstLine;
		[theScanner scanUpToCharactersFromSet: [NSCharacterSet newlineCharacterSet] intoString: &firstLine];
		NSUInteger nDrifters = [[firstLine componentsSeparatedByString:@"\t"] count]-1;
		nDrifters = (nDrifters - 1)/2; // 1 time column, and a column for each x and y
		
        GLDimension *tDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: numberOfLines domainMin: 0 length: numberOfLines];
		GLDimension *drifterDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: nDrifters domainMin: 0 length: nDrifters-1];
        
        self.t = [GLFunction functionOfRealTypeWithDimensions: @[tDim] forEquation:equation];
		self.x = [GLFunction functionOfRealTypeWithDimensions: @[tDim,drifterDim] forEquation:equation];
        self.y = [GLFunction functionOfRealTypeWithDimensions: @[tDim,drifterDim] forEquation:equation];
		NSUInteger i=0;
        theScanner = [NSScanner scannerWithString:fileContents];
        while (theScanner.isAtEnd == NO) {
            NSDecimal val;
            [theScanner scanDecimal:&val];
            NSDecimalNumber *number = [NSDecimalNumber decimalNumberWithDecimal: val];
            self.t.pointerValue[i] = number.doubleValue;
			
			NSUInteger iDrifter = 0;
			while (iDrifter < nDrifters) {
				[theScanner scanDecimal:&val];
				number = [NSDecimalNumber decimalNumberWithDecimal: val];
				self.x.pointerValue[i*nDrifters + iDrifter] = number.doubleValue;
				iDrifter++;
			}
            
            iDrifter = 0;
			while (iDrifter < nDrifters) {
				[theScanner scanDecimal:&val];
				number = [NSDecimalNumber decimalNumberWithDecimal: val];
				self.y.pointerValue[i*nDrifters + iDrifter] = number.doubleValue;
				iDrifter++;
			}
			
            i++;
        }
		

    }
    return self;
}

@end
