#-------------------------------------
# @__{slide()}__
# # Visualization and Visualization in Julia
#
# Turn this file into jupyter notebook with Literate.jl or [jl2nb](https://github.com/j-fu/jl2nb)
#-------------------------------------
# @__{slide()}__
# ##  Plotting & visualization
#
# Human perception is much better adapted to visual representation than to numbers
#
# Purposes of plotting:
# - Visualization of research result for publications & presentations
# - Debugging + developing algorithms
# - "In-situ visualization" of evolving computations
# - Investigation of data
# - 1D, 2D, 3D, 4D data
# - Similar tasks in CAD, Gaming, Virtual Reality $\dots$
# - $\dots$
#-------------------------------------
# @__{slide()}__
# ## Processing steps in visualization
# ### High level tasks:
#
# - Representation of data using elementary primitives: points,lines, triangles, $\dots$
#   - Very different depending on purpose
#
# ### Low level tasks
#
# - Coordinate transformation from "world coordinates" of a particular model to screen coordinates
# - Transformation 3D $\to$ 2D, visibility computation
# - Coloring, lighting, transparency
# - Rasterization: turn smooth data into pixels
#
#-------------------------------------
# @__{slide()}__
# ## Software implementation of low level tasks
#
# - Software: rendering libraries, e.g.  Cairo, AGG
# - Software for vector based graphics formats, e.g. PDF, postscript, svg
# - Typically performed on CPU
#
#-------------------------------------
# @__{slide()}__
# ## Hardware for   low level tasks
# - Huge number of very similar operations 
# - SIMD  parallelism "Single instruction, multiple data" inherent to processing steps in visualization
# - Dedicated hardware: *Graphics Processing Unit* (GPU) frees CPU  from these taks
# - Multiple parallel pipelines, fast memory for intermediate results
#-------------------------------------
# @__{subslide()}__
#
# ### GPU Programming
# - Typically, GPUs are processing units which are connected via bus interface to CPU
# - GPU Programming:
#   - Prepare low level data for GPU
#   - Send data to GPU
#   - Process data in rendering pipeline(s)
# - Modern visualization programs have a CPU part and GPU parts a.k.a. *shaders*
#   - Shaders allow to program details of data processing on GPU
#   - Compiled on CPU, sent along with data to GPU
# - Modern libraries: Vulkan, modern OpenGL/WebGL, DirectX
#-------------------------------------
# @__{fragment()}__
# - Possibility to "mis-use" GPU for numerical computations
# - google "Julia GPU" $\dots$
#-------------------------------------
# @__{subslide()}__
#
# ### GPU Programming in the "old days"
# - "Fixed function pipeline"  in OpenGL 1.1 fixed one particular set of shaders
# - Easy to program
# ````
# glClear()
# glBegin(GL_TRIANGLES)
# glVertex3d(1,2,3)
# glVertex3d(1,5,4)
# glVertex3d(3,9,15)
# glEnd()
# glSwapBuffers()
# ````
# - Not anymore: now you write  shaders for this, compile them, $\dots$
#-------------------------------------
# @__{subslide()}__
# ### Library interfaces to GPU useful for Scientific Visualization
# GPUs are ubiquitous, why aren't they used for visualization ? 
#
# - vtk (backend of Paraview,VisIt)
# - GR framework
# - Three.js (for WebGL in the browser)
# - Makie (as a fresh start in Julia)
# - very few $\dots$
#   - Money seems to be in gaming, battlefield rendering $\dots$
#   - __This is not a Julia only problem__   but also for python, C++, $\dots$
#-------------------------------------
# @__{slide()}__
#
# ## Consequences for Julia
# - It is hard to have high quality and performance for lage datasets at once
# - Julia in many cases (high performance linear algebra for standard float types, sparse matrix solvers $\dots$)
#   relies on well tested libraries
# - Similar approach for graphics
# - [Makie.jl](https://makie.juliaplots.org/stable/): fresh start directly based on OpenGL, but
#   but functionality still behind
#-------------------------------------
# @__{slide()}__
# ### [PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl)
# - Interface to matplotlib from the python world
# - Ability to create publication ready graphs
# - Limited performance (software rendering, many processing steps in python)
# - Best start for users familiar and satisfied with matplotlib performance
#-------------------------------------
# @__{slide()}__
# ### [Plots.jl](https://github.com/JuliaPlots/Plots.jl)
# - Meta package with a number of different backends
# - Backends are already high level libraries 
# - Choose based on  performance, quality, code stability
# - Write code once, just switch backend
#---
# @__{subslide()}__
# 
# #### [GR Framework](https://gr-framework.org/): `Plots.gr()`
# - Design based on [Graphical Kernel System (GKS)](https://en.wikipedia.org/wiki/Graphical_Kernel_System), the first and now nearly forgotten ISO standard for computer graphics as intermediate interface
# - Very flexible concerning low level backend (from  Tektronix to OpenGL...)
# - Corner cases where pyplot has more functionality
# - OpenGL $\Rightarrow$ __fast__! 
# - Few dependencies 
#---
# @__{subslide()}__
# 
# #### [PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl) once again: `Plots.pyplot()`
# - High quality
# - Limited performance
# - Needs python + matplotlib
#
#---
# @__{subslide()}__
# 
# #### More...
#
# - PGFPLots: uses LaTeX typsetting system as backend
#   - probably best quality
#   - slowest (all data are processed via LaTeX)
#   - large dependency 
# - UnicodePlots: ASCII Art revisited
#
#-------------------------------------
# @__{slide()}__
# ## Plots.jl workflow
# - Use fast backend for exploring large data and developing ideas 
# - For creating presentable graphics, prepare data in such a way the they can be quickly loaded 
# - Use high quality backend to tweak graphics for presentation
# - If possible, store graphics in vectorized format
# - See also blog post by [Tamas Papp](https://tamaspapp.eu/post/plot-workflow/)
#----
# @__{subslide()}__
# ### Plots.jl ressources
# - [Learning](http://docs.juliaplots.org/latest/learning/) ressources
# - [Revise.jl](https://github.com/timholy/Revise.jl): automatic file reloading + compilation for REPL based workflow
#-------------------------------------
# @__{slide()}__
# 
# ## Preparing plots
# - We  import `Plots` so you see which methods come from there
import Plots
# - Set a flag variable to check if code runs in jupyter 
injupyter=isdefined(Main, :IJulia) && Main.IJulia.inited
# - A simple plot based on defaults
# - Just use the plot method
function example1(;n=100)
    f(x)=sin(exp(4.0*x))
    X=collect(-1.0:1.0/n:1.0)
    p=Plots.plot(X,f.(X))
end
# Run example
injupyter&&    example1()
#
#-------------------------------------
# @__{slide()}__
# - A good plot has axis labels etc.
# - Also we want to have a better label placement
# - The plot! method allows a successive build-up of the plot using different *attributes*
# - We save the plot to pdf for embedding into presentations
function example2(;n=100)
    f(x)=sin(exp(4x))
    X=collect(-1.0:1.0/n:1.0)
    p=Plots.plot(framestyle=:full,legend=:topleft, title="Example2")
    Plots.plot!(p, xlabel="x",ylabel="y")
    Plots.plot!(p, X,f.(X), label="y=sin(exp(4x))")
    Plots.savefig(p,"example2.pdf")
    return p
end
# Run this example
injupyter&&    example2()
#-------------------------------------
# @__{slide()}__
# - Two functions in one plot
function example3(;n=100)
    f(x)=sin(exp(4x))
    g(x)=sin(exp(-4x))
    X=collect(-1.0:1.0/n:1.0)
    p=Plots.plot(framestyle=:full,legend=:topleft, title="Example3")
    Plots.plot!(p, xlabel="x",ylabel="y")
    Plots.plot!(p, X,f.(X), label="y=sin(exp(4x))", color=Plots.RGB(1,0,0))
    Plots.plot!(p, X,g.(X), label="y=sin(exp(-4x))", color=Plots.RGB(0,0,1))
end
# Run this example
injupyter&&    example3()
#-------------------------------------
# @__{slide()}__
# - Two plots arranged
function example4(;n=100)
    f(x)=sin(exp(4x))
    g(x)=sin(exp(-4x))
    X=collect(-1.0:1.0/n:1.0)
    p1=Plots.plot(framestyle=:full,legend=:topleft, title="Example4")
    Plots.plot!(p1, xlabel="x",ylabel="y")
    Plots.plot!(p1, X,f.(X), label="y=sin(exp(4x))", color=Plots.RGB(1,0,0))
    p2=Plots.plot(framestyle=:full,legend=:topright)
    Plots.plot!(p2, xlabel="x",ylabel="y")
    Plots.plot!(p2, X,g.(X), label="y=sin(exp(-4x))", color=Plots.RGB(0,0,1))
    p=Plots.plot(p1,p2,layout=(2,1))
end
# Run this example
injupyter&&    example4()
#-------------------------------------
# @__{slide()}__
# - Two plots arranged, one scattered
function example5(;n=100)
    f(x)=sin(exp(4x))
    g(x)=sin(exp(-4x))
    X=collect(-1.0:1.0/n:1.0)
    p1=Plots.plot(framestyle=:full,legend=:topleft, title="Example5")
    Plots.plot!(p1, xlabel="x",ylabel="y")
    Plots.plot!(p1, X,f.(X), label="y=sin(exp(4x))", color=Plots.RGB(1,0,0))
    p2=Plots.plot(framestyle=:full,legend=:topright)
    Plots.plot!(p2, xlabel="x",ylabel="y")
    Plots.plot!(p2, X,g.(X), label="y=sin(exp(-4x))", seriestype=:scatter, color=Plots.RGB(0,0,1), markersize=0.5)
    p=Plots.plot(p1,p2,layout=(2,1))
end
# Run this example
injupyter&&    example5()

#-------------------------------------
# @__{slide()}__
# ## Plots terminology
# - Plot: The whole figure/window
# - Subplot: One subplot, containing a title, axes, colorbar, legend, and plot area.
# - Axis: One axis of a subplot, containing axis guide (label), tick labels, and tick marks.
# - Plot Area: The part of a subplot where the data is shown... contains the series, grid lines, etc.
# - Series: One distinct visualization of data. (For example: a line or a set of markers)
# ### Appearance
# - Appearance of Plot, Subplot, Axis, Series is influenced by *attributes*
# - Attributes given as Keyword argument to `Plots.plot()`,`Plots.plot!()`
#-------------------------------------
# @__{subslide()}__
# ## Which attributes are supported in Plots ?
# - "Google" vs. "google the right thing"
#
Plots.plotattr()
#----
# @__{fragment()}__
# $~$
Plots.plotattr(:Series)
#----
# @__{fragment()}__
# $~$
Plots.plotattr("seriestype")
#-------------------------------------
# @__{slide()}__
# - Heatmap with contourlines
function example6(;n=100)
    f(x,y)=sin(10x)*cos(10y)*exp(x*y)
    X=collect(-1.0:1.0/n:1.0)
    Y=view(X,:)
    Z=[f(X[i],Y[j]) for i=1:length(X),j=1:length(Y)]
    p=Plots.plot(X,Y,Z, seriestype=:heatmap,seriescolor=Plots.cgrad([:red,:yellow,:blue]))
    p=Plots.plot!(p,X,Y,Z, seriestype=:contour, seriescolor=:black)
end
# Run this example
injupyter&&    example6()
#-------------------------------------
# @__{slide()}__
# - Heatmap with contourlines plotted during loop
import PyPlot
# $~$
function example7(;n=100,tend=10)
    f(x,y,t)=sin((10+sin(t))*x-t)*cos(10y-t)
    X=collect(-1.0:1.0/n:1.0)
    Y=view(X,:)
    Z=[f(X[i],Y[j],0) for i=1:length(X),j=1:length(Y)]
    for t=1:0.1:tend
        injupyter && IJulia.clear_output(true)
        for i=1:length(X),j=1:length(Y)
            Z[i,j]=f(X[i],Y[j],t)
        end
        p=Plots.plot(title="t=$(t)")
        p=Plots.plot!(p,X,Y,Z, seriestype=:heatmap,seriescolor=Plots.cgrad([:red,:yellow,:blue]))
        p=Plots.plot!(p,X,Y,Z, seriestype=:contour, seriescolor=:black)
        Plots.display(p)
        if Plots.backend_name()==:pyplot
            PyPlot.pause(1.0e-10)
        end
    end
end
# Run this example
injupyter&&    example7()

#-------------------------------------
# @__{slide()}__
# ## Makie
using Makie
# - Create the scene before running
# - Update scene data via lift in the inner loop
function example8(;n=100,tend=10)
    scene=Makie.Scene()
    f(x,y,t)=sin((10+sin(t))*x-t)*cos(10y-t)
    X=collect(-1.0:1.0/n:1.0)
    Y=view(X,:)
    Z=[f(X[i],Y[j],0) for i=1:length(X),j=1:length(Y)]
    scene_data=Makie.Node(Z)
    scene_title=Makie.Node("t=0.0")
    scene_xcoord=Makie.Node(X)
    scene_ycoord=Makie.Node(X)
    Makie.heatmap!(scene, X,Y,lift(a->a,scene_data),levels=10)
    Makie.contour!(scene, X,Y,lift(a->a,scene_data),levels=10,color=:black)
    Makie.display(scene)

    for t=1:0.1:tend
        for i=1:length(X),j=1:length(Y)
            Z[i,j]=f(X[i],Y[j],t)
        end
        push!(scene_data,Z)
        Makie.sleep(1.0e-10)
    end
end
# Run this example
injupyter&&    example8()


#-------------------------------------
# @__{slide()}__
# ## [VTK](http://vtk.org)
# - Visualization primitives in scientific computing:
#   - Datasets on rectangular and unstructured discretization grids
#   - Scalar data
#   - Vector data
# - *Visualization Toolkit* vtk provides  an API based on these
#     primitives and uses up-to date graphics API (OpenGL) to render these data
#   - Well maintained, "working horse" in  high performance computing
#   - Open Source
#   - Paraview, VisIt: GUI programs around vtk
#   - APIs for C++, Python, $\dots$
#   - Program the piplines you usually create via GUI in paraview
# - [WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl): Julia package for writing VTK data files
#-------------------------------------
# @__{subslide()}__
# ## VTKFig
# 
# [VTK](http://vtk.org)  based C++ graphics library for plotting and for data on rectilinear and unstructured grids with an flexible and easy to use API.
# 
# ## Features
# 
# - Standard views for scalars, vectors and grids in 2D, 3D:
#   - isolines, filled contours
#   - plane cuts, isosurfaces
#   - quiver plots, stream ribbons
# - Separate  rendering thread  allowing for  handling of  changing data   managed by the vtkfig::Frame class
# - Extensible by implementing derived classes  containing  vtk rendering pipelines from vtkfig::Figure
# - (Emerging) C API
# - Experimental companion Julia package [VTKFig.jl](https://github.com/j-fu/VTKFig.jl) based on C API and Julia `ccall`
#-------------------------------------
# @__{slide()}__
# - Add unregistered VTKFig.jl ( `Pkg.add("https://github.com/j-fu/VTKFig.jl")` )
# - Checkout and compile vtkfig
# - Have vtkfig.so on LD_LIBRARY_PATH
# - Heatmap with contourlines plotted during loop based on vtkfig
import VTKFig
# $~$
function example9(;n=100,tend=10)
    f(x,y,t)=sin((10+sin(t))*x-t)*cos(10y-t)
    X=collect(-1.0:1.0/n:1.0)
    V=[f(X[i],X[j],0) for i=1:length(X),j=1:length(X)]
    VTKFig.clear()
    griddata=VTKFig.DataSet()
    contour=VTKFig.ScalarView()
    VTKFig.set_size(500,500)
    VTKFig.set_rectilinear_grid(griddata,X,X) 
    VTKFig.set_point_scalar(griddata,vec(V),"V")
    VTKFig.set_data(contour,griddata,"V")
    VTKFig.add_figure(contour)
    for t=1:0.1:tend
        for i=1:length(X),j=1:length(X)
            V[i,j]=f(X[i],X[j],t)
        end
        VTKFig.set_frame_title("t=$(t)")
        VTKFig.set_point_scalar(griddata,vec(V),"V")
        VTKFig.show()
    end
end
# Run this example
injupyter&&    example9()

#-------------------------------------
# @__{slide()}__
# - Heatmap with contourlines plotted during loop based on vtkfig
function example10(;n=50,tend=10)
    f(x,y,z,t)=sin((10+sin(t))*x-t)*cos(10y-t)*sin(3z-2t)
    X=collect(-1.0:1.0/n:1.0)
    N=length(X)
    V=[f(X[i],X[j],X[k],0) for i=1:N,j=1:N,k=1:N]
    VTKFig.clear()
    griddata=VTKFig.DataSet()
    contour=VTKFig.ScalarView()
    VTKFig.set_size(500,500)
    VTKFig.set_rectilinear_grid(griddata,X,X,X) 
    VTKFig.set_point_scalar(griddata,vec(V),"V")
    VTKFig.set_data(contour,griddata,"V")
    VTKFig.show_isosurfaces(contour,true)
    VTKFig.set_isolevels(contour,[0.1,0.9])
    VTKFig.add_figure(contour)
    for t=1:0.1:tend
        for i=1:N,j=1:N,k=1:N
            V[i,j,k]=f(X[i],X[j],X[k],t)
        end
        VTKFig.set_frame_title("t=$(t)")
        VTKFig.set_point_scalar(griddata,vec(V),"V")
        VTKFig.show()
    end
end
# Run this example
injupyter&&    example10()

