using Documenter
using IntervalCensoredMultistate

makedocs(
  sitename="IntervalCensoredMultistate.jl",
  modules=[IntervalCensoredMultistate],
  format=Documenter.HTML(
    prettyurls=get(ENV, "CI", "false") == "true",
    canonical="https://luyouepiusf.github.io/IntervalCensoredMultistate.jl/"
  ),
  pages=[
    "Home" => "index.md",
    "API Reference" => "api.md",
  ],
  checkdocs=:none,
)

deploydocs(
  repo="github.com/luyouepiusf/IntervalCensoredMultistate.jl.git",
  devbranch="master",
)
