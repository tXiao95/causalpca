-- Remove rendered DOI links without changing the shared bibliography data.
local function strip_doi(div)
  if div.identifier ~= "refs" then
    return nil
  end
  return div:walk({
    Inlines = function(inlines)
      local result = pandoc.List()
      local after_doi = false
      for _, inline in ipairs(inlines) do
        if inline.t == "Link" and inline.target:match("doi%.org/") then
          if #result > 0 and result[#result].t == "Space" then
            result:remove(#result)
          end
          after_doi = true
        elseif after_doi and inline.t == "Str" and inline.text == "." then
          after_doi = false
        else
          result:insert(inline)
          after_doi = false
        end
      end
      return result
    end,
  })
end

function Pandoc(doc)
  return pandoc.utils.citeproc(doc):walk({ Div = strip_doi })
end
