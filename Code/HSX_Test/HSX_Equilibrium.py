import vmecpp


'''With Bootstrap Currents'''
model_file = '../../Data/FLARE_DB/HSX_Test/shared_data/input.hsx'
vmec_input = vmecpp.VmecInput.from_file(model_file) 

vmec_output = vmecpp.run(vmec_input)


print(vmec_output.mercier.iota)
vmec_output.wout.save("wout_w7x_curr.nc")