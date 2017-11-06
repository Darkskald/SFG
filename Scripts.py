from ipy_interpreter import Ipy_Interpreter as ip
from Classes import Plotter

I = ip()

sens = ["BX9", "BX12"]


for sensitizer in sens:

    s = "su "+sensitizer
    I.get(s)
    dates = []
    for spectrum in I.subset:
        date = spectrum.name.date
        if date not in dates:
            dates.append(date)

        backup = I.subset
        for date in dates:

            I.refine("d "+date)
            samples = []

            for spectrumx in I.subset:
                sample = spectrumx.name.sample_number
                if sample not in samples:
                    samples.append(sample)

            for sample in samples:
                I.refine("s "+str(sample))
                t = sensitizer+" "+str(date)+"_"+str(sample)
                P = Plotter(I.subset,title=t)
                P.custom_plot()
                I.recovery()

            I.subset = backup


