using Documenter
using ARMtools

makedocs(
    sitename = "ARMtools",
    format = Documenter.HTML(),
    modules = [ARMtools]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
