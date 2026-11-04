-- Use the slide theme's dark blue for inline citations in PDF output.
function Cite(citation)
  if FORMAT:match("latex") or FORMAT:match("beamer") then
    return {
      pandoc.RawInline("latex", "{\\color{emoryblue}"),
      citation,
      pandoc.RawInline("latex", "}"),
    }
  end
end
