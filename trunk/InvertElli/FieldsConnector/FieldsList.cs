using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace FieldsConnector
{
    public class FieldsList
    {

        private Dictionary<string, AbstractFieldConnector> _list;
        public AbstractFieldConnector this[string index]
        {

            get
            {
                return _list[index];

            }

        }
        public void addTier(AbstractFieldConnector afc)
        {
            _list.Add(afc.Name, afc);
        }

    }
}
