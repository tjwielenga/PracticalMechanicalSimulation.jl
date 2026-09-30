-- Number displayed equations and resolve symbolic equation references.
--
-- An unlabeled displayed equation needs no special markup:
--
--   $$
--   a=b
--   $$
--
-- Add a stable label only when the equation must be referenced:
--
--   $$
--   a=b
--   $$ {#eq:example}
--
-- Refer to the label as `Eq. @eq:example` or `Equation @eq:example`.

local equation_numbers = {}
local equation_count = 0
local output_equation_count = 0

local function split_display_math(para)
    if para.tag ~= "Para" then
        return nil
    end

    local has_display_math = false
    for _, inline in ipairs(para.content) do
        if inline.tag == "Math" and inline.mathtype == "DisplayMath" then
            has_display_math = true
            break
        end
    end
    if not has_display_math then
        return nil
    end

    local blocks = {}
    local ordinary = pandoc.Inlines({})
    local function flush_ordinary()
        while #ordinary > 0 and
                (ordinary[#ordinary].tag == "Space" or
                 ordinary[#ordinary].tag == "SoftBreak") do
            ordinary:remove(#ordinary)
        end
        if #ordinary > 0 then
            table.insert(blocks, pandoc.Para(ordinary))
            ordinary = pandoc.Inlines({})
        end
    end

    local index = 1
    while index <= #para.content do
        local inline = para.content[index]
        if inline.tag == "Math" and inline.mathtype == "DisplayMath" then
            flush_ordinary()
            local equation = pandoc.Inlines({inline})
            local marker_index = index + 1
            if para.content[marker_index] and
                    para.content[marker_index].tag == "Space" then
                equation:insert(para.content[marker_index])
                marker_index = marker_index + 1
            end
            local marker = para.content[marker_index]
            if marker and marker.tag == "Str" and
                    marker.text:match("^%{#eq:[%w_.:-]+%}$") then
                equation:insert(marker)
                index = marker_index
            end
            table.insert(blocks, pandoc.Para(equation))
        elseif inline.tag ~= "SoftBreak" or #ordinary > 0 then
            ordinary:insert(inline)
        end
        index = index + 1
    end
    flush_ordinary()
    return blocks
end

local function displayed_equation(para)
    if para.tag ~= "Para" or #para.content == 0 then
        return nil, nil
    end

    local math = para.content[1]
    if math.tag ~= "Math" or math.mathtype ~= "DisplayMath" then
        return nil, nil
    end

    if #para.content == 1 then
        return math, nil
    end

    local marker_index = 2
    if para.content[2].tag == "Space" then
        marker_index = 3
    end
    if marker_index ~= #para.content or
            para.content[marker_index].tag ~= "Str" then
        return nil, nil
    end

    local label = para.content[marker_index].text:match("^%{#(eq:[%w_.:-]+)%}$")
    if label == nil then
        return nil, nil
    end
    return math, label
end

local function register_equation(para)
    local math, label = displayed_equation(para)
    if math == nil then
        return
    end

    equation_count = equation_count + 1
    if label ~= nil then
        if equation_numbers[label] ~= nil then
            error("duplicate equation label '" .. label .. "'")
        end
        equation_numbers[label] = equation_count
    end
end

local function equation_table(math, number, label)
    local empty_cell = pandoc.Cell({}, pandoc.AlignLeft, 1, 1)
    local equation_cell = pandoc.Cell(
        {pandoc.Plain({math})}, pandoc.AlignCenter, 1, 1)
    local number_cell = pandoc.Cell(
        {pandoc.Plain({pandoc.Str("(" .. number .. ")")})},
        pandoc.AlignRight, 1, 1)
    local row = pandoc.Row({empty_cell, equation_cell, number_cell})
    local body = pandoc.TableBody({row}, {}, 0)
    local identifier = label or ("equation-" .. number)

    return pandoc.Table(
        pandoc.Caption(),
        {
            {pandoc.AlignLeft, 0.08},
            {pandoc.AlignCenter, 0.84},
            {pandoc.AlignRight, 0.08},
        },
        pandoc.TableHead({}),
        {body},
        pandoc.TableFoot({}),
        pandoc.Attr(identifier, {"numbered-equation"})
    )
end

local function numbered_equation(para)
    local math, label = displayed_equation(para)
    if math == nil then
        return nil
    end

    local number
    if label ~= nil then
        number = equation_numbers[label]
    else
        number = nil
    end

    if number == nil then
        output_equation_count = output_equation_count + 1
        number = output_equation_count
    else
        output_equation_count = output_equation_count + 1
        if number ~= output_equation_count then
            error("internal equation-numbering order mismatch")
        end
    end

    if FORMAT:match("latex") then
        local latex = "\\begin{equation}\\tag{" .. number .. "}"
        if label ~= nil then
            latex = latex .. "\\label{" .. label .. "}"
        end
        latex = latex .. "\n" .. math.text .. "\n\\end{equation}"
        return pandoc.RawBlock("latex", latex)
    end

    return equation_table(math, number, label)
end

local function equation_reference(cite)
    if #cite.citations ~= 1 then
        return nil
    end

    local label = cite.citations[1].id
    if not label:match("^eq:") then
        return nil
    end

    local number = equation_numbers[label]
    if number == nil then
        error("undefined equation label '" .. label .. "'")
    end

    return pandoc.Link(
        {pandoc.Str("(" .. number .. ")")},
        "#" .. label,
        "Equation " .. number
    )
end

function Pandoc(document)
    equation_numbers = {}
    equation_count = 0
    output_equation_count = 0

    document = document:walk({Para = split_display_math})
    document:walk({Para = register_equation})
    return document:walk({
        Cite = equation_reference,
        Para = numbered_equation,
    })
end
