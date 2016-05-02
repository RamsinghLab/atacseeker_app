{% comment %}
    
    Welcome to the BaseSpace Report Builder!

{% endcomment %}

<td>{{ result.files["atacseeker.html"] | append: "bar" }}</td>

{% for key in result.files %}
<iframe src={{ result.files[key].href }} frameborder="0" style="overflow:hidden;height:100%;width:100%" height="100%" width="100%">
</iframe>
{% endfor %}
