from browse.search import search_types
import json
from django import template
register = template.Library()

@register.filter('fieldtype')
def fieldtype(field):
    return field.field.__class__.__name__

@register.filter('bootchoice_choice')
def bootchoice_choice(field):
    #choice = str(field).replace("<select", "<select width=100% class=\"form-control\"")
    #return choice

    if field.label != "Histone":
        html = get_search_type(field)
        html += '<input type="text" class="form-control" name="{0}" id="{0}" />'.format(field.id_for_label)
        html += get_pull_down(field.field.queryset.values_list("id", flat=True), field.id_for_label, reset="text")
    else:
        html = get_pull_down(field.field.queryset.values_list("id", flat=True), field.id_for_label, default_name="----")

    return html

@register.filter('simple_choice')
def simple_choice(field):
    html = '<input type="hidden" name="{0}" id="{0}" value="{1}">'.format(field.id_for_label, field.field.initial)
    html += '<div class="btn-group" role="group" aria-label="..." id="{}_button_group">'.format(field.id_for_label)
    for choice in field.field.choices:
        color = "primary" if field.field.initial == choice[0] else "default"
        html += '<button type="button" class="btn btn-{} data-form-{}" data-value="{}">{}</button>'.format(color, field.id_for_label, choice[0], choice[1])
    html += '</div>'
    js = '<script type="text/javascript">'
    js += '$(".data-form-{}").on("click", function() {{'.format(field.id_for_label)
    js +=     '$("#{}_button_group").find("button").removeClass("btn-primary").addClass("btn-default");'.format(field.id_for_label)
    js +=     '$(this).removeClass("btn-default").addClass("btn-primary");'
    js +=     '$("#{}").val($(this).attr("data-value"));'.format(field.id_for_label)
    js += "});"
    js += '</script>' 
    html += js
    return html

@register.filter('get_search_type')
def get_search_type(field):
    if fieldtype(field) in ["CharField", "ModelChoiceField"]:
        search_type = search_types[str]
    elif fieldtype(field)in ["IntegerField", "FloatField"]:
        search_type = search_types[int]
    else:
        return ""
    return get_pull_down(search_type.keys(), field.id_for_label+"_search_type")

def get_pull_down(names, id, reset="menu", default_name=""):
    """This will turn a list into a pulldown menu from Bootstrap, Assoociated javascript is also
    created to change the form value

    NOTE: This is extrmeley hacky and will be replaced by django-bootsrap3
    """
    html = '<div class="input-group-btn">'
    if reset == "menu":
        default_name = default_name or "is"
        html += '<input type="hidden" name="{0}" id="{0}" value="{1}">'.format(id, default_name if default_name == "is" else "")
    else:
        default_name= ""
    html += '<button type="button" id="{}_drop_down_button" class="btn btn-default dropdown-toggle" data-toggle="dropdown" aria-expanded="false" width="100%">{} <span class="caret"></span></button>'.format(id, default_name)
    html += '<ul class="dropdown-menu" role="menu" id="{}_drop_down">'.format(id)
    js = '<script type="text/javascript">'
    if default_name != "is":
        html +=  '<li><a href="#" id="{}_drop_down_default">{}</a></li>\n'.format(id, default_name)
        js += "$('#{}_drop_down_default').on('click', function(){{ ".format(id)
        js += "$('#{}_drop_down_button').html($(this).html() + ' <span class=\"caret\"></span>'); ".format(id)
        js +=  "$('#{}').val(""); ".format(id)
        js += "});"
    for i, name in enumerate(names):
        html +=  '<li><a href="#" id="{}_drop_down_{}">{}</a></li>\n'.format(id, i, name)
        js += "$('#{}_drop_down_{}').on('click', function(){{ ".format(id, i)
        if reset == "menu":
            js +=     "$('#{}_drop_down_button').html($(this).html() + ' <span class=\"caret\"></span>'); ".format(id)
        js +=         "$('#{}').val($(this).html()); ".format(id)
        js += "});"
    js += "</script>"
    html += '</ul></div>'
    html += js
    return html

@register.filter('jsonify')
def jsonify(object):
    return json.dumps(object)

@register.filter('bootstrapify')
def bootstrapify(field):
    field_types = {"CharField":"text", "IntegerField":"int", "ModelForm":None}
    field_type = field.field.__class__.__name__
    try:
        search_type = field_types[field_type]
    except KeyError:
        return ""

    html += '<div class="form-group">'
    html += field.label_tag
    html += '<div class="input-group">'
    html +=    '<div class="input-group-btn">'
    if search_type in search_types:
        #js = '<script type="text/javascript">'
        #html +=    '<button type="button" class="btn btn-default dropdown-toggle" data-toggle="dropdown" aria-expanded="false">Is <span class="caret"></span></button>'
        #html +=    '<input type="hidden" name="{0}_search_type" id="{0}_search_type" value="is">'.format(field.id_for_label)
        #html +=    '<ul class="dropdown-menu" role="menu" id="{}_search_type_drop_down">'.format(field.id_for_label)
        #for i, name in enumerate(search_types[search_type].keys()):
        #    html +=  '<li><a href="#" id="{}_{}_search_type_drop_down">{}</a></li>\n'.format(field.id_for_label, i, name)
        #    js += "$('#{}_{}_search_type_drop_down').on('click', function(){ ".format(field.id_for_label, i)
        #    js +=     "$('#{}_search_type_drop_down').html($(this).html() + '<span class=\"caret\"></span>'); "
        #    js +=     "$('#{}_search_type').val($(this).html()); "
        #    js += "});"
        #html +=    '</ul>'
        #js += "</script>"
        #jtml += js
        html +=    '<input type="{}" class="form-control" id="{}" />'.format(field_types[search_type], field.id_for_label)
        html +='</div><!-- /btn-group -->'
    elif field_type == "ModelForm":
        html += str(field).replace("<select", "<select class=\"form-control\"")
    if field.help_text:
        html +='<span class="help-text">{}</span>'.format(field.help_text)
    html += '</div>'

    return html