-- filtres/titre-environnements.lua
local compteurs = {}

local etiquettes = {
  ["env-definition"] = "Définition",
  ["env-propriete"]  = "Propriété",
  ["env-remarque"]   = "Remarque",
  ["env-exemple"]    = "Exemple",
  ["env-theo"]    = "Pour aller plus loin avec la théorie"
}


function Div(el)
  for classe, etiquette in pairs(etiquettes) do
    if el.classes:includes(classe) then
      local texte = etiquette
      if el.attributes["titre"] then
        texte = texte .. " — " .. el.attributes["titre"]
      end

      local inlines = pandoc.read(texte, "markdown").blocks[1].content

      local titre = pandoc.Para({
        pandoc.Span(inlines, pandoc.Attr("", {"env-titre"}))
      })

      table.insert(el.content, 1, titre)
      return el
    end
  end
end