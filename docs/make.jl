using Documenter
using CFSjul

makedocs(
    sitename = "CFSjul Documentation",
    pages = [
        "Index" => "index.md",
        "An other page" => "anotherPage.md",
    ],
    format = Documenter.HTML(),
    modules = [CFSjul]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "git@github.com:iperezav/CFSjul.git",
    devbranch = "main"
)
