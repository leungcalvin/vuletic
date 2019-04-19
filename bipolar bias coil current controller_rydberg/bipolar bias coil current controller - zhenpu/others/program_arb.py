import visa
import numpy as np
import serial

rm = visa.ResourceManager('@py')
visa_addr = "Agilent_Arb"
print "Opening connection to VISA address %s" % visa_addr
inst = rm.open_resource(visa_addr)

            
    
def DC(V):
    """ set output to DC mode at voltage V """
    
    inst.write("APPL:DC DEF,DEF, %.4f" % V)
    inst.write("OUTP ON")
        
def _sendWf(wf_data,Vbot,Vtop,ttot,ncycles=1,load=50):
        # function for internal use, takes array of data wf_data,
        # min,max voltages (in Volts)
        # total time of sequence (in seconds)
        # number of cycles
        # and impedance of output (to interpret voltages)
        
        # wf_data should range from -2047 to 2047
        
        #print wf_data

        #nBytesStr = str(2*array_len)
        #byteStrLen = len(nBytesStr)

        command_string = ":OUTP OFF;:DATA:DAC VOLATILE"
        for pt in wf_data:
            command_string += ", %d" % int(pt)
            
        inst.write(command_string)

      
        
        command_strings = []
        command_strings.append(r':FUNC:USER VOLATILE')
        command_strings.append(r':FUNC:SHAP USER')
        command_strings.append(r':TRIG:SOUR IMM')
        command_strings.append(r':BURS:STAT ON')
        command_strings.append(r':BURS:NCYC '+str(ncycles))
        command_strings.append(r':BURS:MODE TRIG')
        command_strings.append(r':BURS:INT:PER '+str(ttot+0.00000001))

        # important to set this before specifying voltages
        if load>50:
            command_strings.append(r':OUTP:LOAD INF')
        else:
            command_strings.append(r':OUTP:LOAD INF')
        command_strings.append(r':OUTP:POL NORM')

    
        command_strings.append(r':VOLT:RANG:AUTO ON')

        freq = 1.0/ttot
        command_strings.append(':FREQ %f' % freq)

        command_strings.append(r':VOLT:HIGH %f' % Vtop)
        command_strings.append(r':VOLT:LOW %f' % Vbot)

        command_strings.append(r':OUTP ON')

        command_total = ''
        for s in command_strings:
            command_total += (s + ';')
            
        inst.write(command_total)

        print "Sent: %s" % command_total
        

def arb(total_time, voltages):
        #The voltage values should go from -2047 to +2047.
        # the Vmin/Vmax commands specify how these map to actual output voltages.    
        voltage_max = max(voltages)
        voltage_min = min(voltages)
        data_scaled = [ int(round(2*2047*(x - voltage_min)/(voltage_max - voltage_min) - 2047)) for x in voltages]
        
        _sendWf(data_scaled,voltage_min,voltage_max,total_time)
        
        return 1        

x = np.linspace(0,1,11)
y = np.sin(10*x)
arb(0.003,x)