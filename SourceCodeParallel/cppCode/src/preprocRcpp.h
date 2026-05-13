#ifndef PREPROC_H
#define PREPROC_H

// parameters to define prior to compile time

#define UseParallel 					@ADDuseParallel
#define UseEventTime 					@ADDuseEventTime
#define UseEventVar 					@ADDuseEventVar
#define ReturnRK4 						@ADDreturnRK4
#define VerboseCMAES			        @ADDVerboseCMAES

#define ntForPython						@ADDNtForPython
#define TInitForPython					@ADDTInitForPython
#define TEndForPython					@ADDTEndForPython
#define NVForPython                     @ADDNVForPython
#define NIVForPython                    @ADDNIVForPython
#define VarNamesForPython				@ADDVarNamesForPython
#define IntermediateVarNamesForPython	@ADDIntermediateVarNamesForPython
#define ParmsNamesForPython				@ADDParmsNamesForPython
#define YInitForPython					@ADDYInitForPython
#define ParmsForPython					@ADDParmsForPython
#define SamplesExogVarForPython			@ADDSamplesExogVarForPython
#define NSamplesVarExogVarForPython		@ADDNSamplesVarExogVarForPython
#define NVarExogVarForPython			@ADDNVarExogVarForPython


typedef @ADDTInt TInt;
#endif
