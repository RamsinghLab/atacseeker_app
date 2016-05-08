{% comment %}
    
    Welcome to the BaseSpace Report Builder!

{% endcomment %}


{% for key in result.files %}
    {% if key contains 'atacseeker.html' %}
        <iframe src={{ result.files[key].href }} frameborder="0" style="overflow:hidden;height:100%;width:100%" height="100%" width="100%">
        </iframe>
    {% endif %}
{% endfor %}