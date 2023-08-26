using PlotlyJS
using Colors


function set_figure_style!(p; title = "", ylabel = "", xlabel ="")
    relayout!(p,
        title=title, 
        xaxis=attr(
            title=xlabel, 
            showline=true,
            linecolor="black",
            mirror=true,
            ticks="outside",
            gridcolor="lightgrey"
        ),
        yaxis=attr(
            title=ylabel, 
            showline=true,
            linecolor="black",
            mirror=true,
            ticks="outside",
            gridcolor="lightgrey"
        ),
        plot_bgcolor="white",
        width=1000,
        height=400
    )
end
function plotZeemanMap(mol::MoleculeHamiltonian, sol_vec; N = [0], energyRef = 0)
    basisUC = getBasisUC(mol.MolOp.basisTree)
    spinDim = prod([length(NuclearSpin.spin) for NuclearSpin in endNode(mol.MolOp.basisTree)[2:end]])
    indsOI = [[sum([(2*N_i2 + 1)*spinDim for N_i2 in 0:(N_i-1)]), sum([(2*N_i2 + 1)*spinDim for N_i2 in 0:N_i])] for N_i in N]
    p = plot()
    for rotationalStates in indsOI,  state in (rotationalStates[1] + 1):rotationalStates[2]
        addtraces!(p, scatter(x = [sol_i.B_field*1e4 for sol_i in sol_vec], y = [sol_i.val[state] for sol_i in sol_vec]*1e-9 .- energyRef*1e-9, text=[KetName(sol_i.vec[:, state], basisUC) for sol_i in sol_vec], name = "State $state"))
    end
    set_figure_style!(p, title = "Zeeman Scan", xlabel = "B-Field (G)", ylabel = "Energy (GHz)")
    p
end

function plotStarkMap(mol::MoleculeHamiltonian, StarkScan; N = [0], energyRef = 0)
    basisUC = getBasisUC(mol.MolOp.basisTree)
    spinDim = prod([length(NuclearSpin.spin) for NuclearSpin in endNode(mol.MolOp.basisTree)[2:end]])
    indsOI = [[sum([(2*N_i2 + 1)*36 for N_i2 in 0:(N_i-1)]), sum([(2*N_i2 + 1)*36 for N_i2 in 0:N_i])] for N_i in N]
    p = plot()
    for rotationalStates in indsOI,  state in (rotationalStates[1] + 1):rotationalStates[2]
        
        addtraces!(p, scatter(x = [sol_i.E_field*1e-2 for sol_i in StarkScan], y = [sol_i.val[state] for sol_i in StarkScan]*1e-9 .- energyRef*1e-9,  text=[KetName(sol_i.vec[:, state], basisUC) for sol_i in StarkScan], name = "State $state"))
    end
    set_figure_style!(p, title = "Stark Scan", xlabel = "E-Field (V/cm)", ylabel = "Energy (GHz)")
    p
end


function plotIntensityScan(mol::MoleculeHamiltonian, sol_vec; N = [0], Beam = 1, energyRef = 0)
    basisUC = getBasisUC(mol.MolOp.basisTree)
    spinDim = prod([length(NuclearSpin.spin) for NuclearSpin in endNode(mol.MolOp.basisTree)[2:end]])
    indsOI = [[sum([(2*N_i2 + 1)*spinDim for N_i2 in 0:(N_i-1)]), sum([(2*N_i2 + 1)*spinDim for N_i2 in 0:N_i])] for N_i in N]
    p = plot()
    for rotationalStates in indsOI,  state in (rotationalStates[1] + 1):rotationalStates[2]
        addtraces!(p, scatter(x = [sol_i.Intensity[Beam] for sol_i in sol_vec], y = [sol_i.val[state]  for sol_i in sol_vec]*1e-9 .- energyRef*1e-9, text=[KetName(sol_i.vec[:, state], basisUC) for sol_i in sol_vec], name = "State $state"))
    end
    set_figure_style!(p, title = "Intensity Scan", xlabel = "Intensity (W/m)", ylabel = "Energy (GHz)")
    p
end
